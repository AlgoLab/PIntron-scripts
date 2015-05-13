#!/usr/bin/env python

from __future__ import print_function

import argparse
import contextlib
import gzip
import json
import logging
import re
import sys


@contextlib.contextmanager
def smart_open_out(filename=None):
    if filename and filename != '-':
        fn = filename
        fo = open(filename, 'w')
        if filename.endswith(".gz"):
            fh = gzip.GzipFile(fn, 'wb', 9, fo)
        else:
            fh = fo
    else:
        fn = "<stdout>"
        fo = sys.stdout
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not fo:
            fh.close()
        if fo is not sys.stdout:
            fo.close()


def parse_est_fasta_header(line):
    identifier = line
    if '/gb=' in identifier:
        identifier = re.search('\/gb=([A-Z_0-9:]+)', identifier).groups()[0]
    return identifier


def parse_block_alignment_str(line, genomic_block):
    parts = line.split()
    block = {'transcript_start': int(parts[0]),
             'transcript_end': int(parts[1]),
             'genomic_relative_start': int(parts[2]),
             'genomic_relative_end': int(parts[3]),
             'transcript_sequence': parts[4],
             'genomic_sequence': parts[5]}
    if genomic_block['strand'] == '+':
        block['genomic_absolute_start'] = genomic_block['start'] + block['genomic_relative_start']
        block['genomic_absolute_end'] = genomic_block['start'] + block['genomic_relative_end']
    else:
        block['genomic_absolute_start'] = genomic_block['end'] - block['genomic_relative_start'] + 1
        block['genomic_absolute_end'] = genomic_block['end'] - block['genomic_relative_end'] + 1
    return block


def pintron_alignments(alignments, genomic_block):
    curr_alignment = None
    for line in alignments:
        line = line.strip()
        if line.startswith(">"):      # New sequence/alignment
            if curr_alignment:
                yield curr_alignment
            curr_alignment = {'identifier': parse_est_fasta_header(line),
                              'blocks': []}
        elif line.startswith("#"):    # Annotation (currently ignored)
            pass
        elif line[0].isdigit():       # New block
            assert(curr_alignment)
            curr_alignment['blocks'].append(parse_block_alignment_str(line, genomic_block))
    if curr_alignment:
        yield curr_alignment


def parse_genomic_header(line):
    parts = line.split(":")
    strand = int(parts[3]+"1")
    strand = "+" if strand > 0 else "-"
    return {'seqname': parts[0],
            'start': int(parts[1]),
            'end': int(parts[2]),
            'strand': strand}


def parse_introns(introns_stream, genomic_block, alignments):
    introns = []
    for identifier, line in enumerate(introns_stream):
        intron = {"identifier": "Int{0}".format(identifier)}
        (intron['relative_start'], intron['relative_end'],
         intron['absolute_start'], intron['absolute_end'],
         intron['length'],
         intron['number_of_supporting_transcripts'], seq_list,
         intron['donor_alignment_error'], intron['acceptor_alignment_error'],
         intron['donor_score'], intron['acceptor_score'],
         intron['BPS_score'], intron['BPS_position'],
         intron['type'], intron['pattern'], intron['repeat_sequence'],
         intron['donor_exon_suffix'], intron['prefix'], intron['suffix'],
         intron['acceptor_exon_prefix']) = line.rstrip().split("\t")
        intron["supporting_transcripts"] = {i: {} for i in seq_list.split(',') if i != ''}

        for field in ('relative_start', 'relative_end',
                      'absolute_start', 'absolute_end',
                      'length', 'number_of_supporting_transcripts',
                      'BPS_position'):
            intron[field] = int(intron[field])

        # a bug in PIntron (file predicted-introns.txt) causes to report absolute
        # coordinates of introns on strand '+' shifted to the left of 1bp
        if genomic_block['strand'] == "+":
            intron['absolute_start'] = intron['absolute_start'] + 1
            intron['absolute_end'] = intron['absolute_end'] + 1

        for field in ('donor_alignment_error', 'acceptor_alignment_error',
                      'donor_score', 'acceptor_score', 'BPS_score'):
            intron[field] = float(intron[field])

        if intron['BPS_position'] < 0:
            del intron['BPS_position']
        introns.append(intron)

    # for each intron, add the alignment of the surrounding exons.
    # Since different factorizations can support the same intron, the first
    # step is to find all pairs of exons supporting an intron
    def supporting_factors(intron):
        pairs = []
        for identifier in intron["supporting_transcripts"].keys():
            factor = alignments[identifier]
            good_left = [exon for exon in factor['blocks']
                         if exon['genomic_relative_end'] == intron['relative_start'] - 1]
            good_right = [exon for exon in factor['blocks']
                          if exon['genomic_relative_start'] == intron['relative_end'] + 1]
            if len(good_left) == 1 and len(good_right) == 1:
                pairs.append((identifier, good_left[0], good_right[0]))
            else:
                logging.error("Could not find a supporting alignment of " +
                              "transcript %s for intron %s",
                              identifier, intron['identifier'])
        assert len(pairs) == intron['number_of_supporting_transcripts']
        return pairs

    #
    # Each intron has the list of supporting_transcripts.
    # For each such seq we provide the suffix/prefix of the prev/next exon
    for intron in introns:
        for [seq, donor_factor, acceptor_factor] in supporting_factors(intron):
            intron['supporting_transcripts'][seq] = {
                'donor_factor_suffix':
                donor_factor['transcript_sequence'][-len(intron['donor_exon_suffix']):],
                'acceptor_factor_prefix':
                acceptor_factor['transcript_sequence'][:len(intron['acceptor_exon_prefix'])],
                'acceptor_factor_start': acceptor_factor['transcript_start'],
                'donor_factor_end': donor_factor['transcript_end'],
                'acceptor_factor_end': acceptor_factor['transcript_end'],
                'donor_factor_start': donor_factor['transcript_start'],
            }

    return introns


def convert_to_dict(genomic_file, alignment_file, introns_file):
    logging.debug("Parsing genomic coordinates from file '%s'", genomic_file.name)
    genomic_block = parse_genomic_header(genomic_file.readline().strip().lstrip(">"))
    logging.debug("Read: %s", genomic_block)

    logging.debug("Parsing alignments from file '%s'", alignment_file.name)
    alignments = {alignment['identifier']: alignment
                  for alignment in pintron_alignments(alignment_file, genomic_block)}
    logging.debug("Read %d alignments", len(alignments))

    logging.debug("Parsing predicted introns from file '%s'", introns_file.name)
    introns = parse_introns(introns_file, genomic_block, alignments)
    logging.debug("Read %d introns", len(introns))

    return {'genome': genomic_block,
            'introns': introns,
            'alignments': alignments}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--pintron-genomic-file',
                        nargs='?',
                        type=argparse.FileType(mode='r'),
                        default=open('genomic.txt'))
    parser.add_argument('-a', '--pintron-align-file',
                        nargs='?',
                        type=argparse.FileType(mode='r'),
                        default=open('out-after-intron-agree.txt'))
    parser.add_argument('-i', '--pintron-introns-file',
                        nargs='?',
                        type=argparse.FileType(mode='r'),
                        default=open('predicted-introns.txt'))
    parser.add_argument('-j', '--output-json-file',
                        nargs='?',
                        default="-")
    parser.add_argument('-v', '--verbose',
                        help='increase output verbosity',
                        action='count', default=0)

    args = parser.parse_args()
    args = vars(args)
    if args['verbose'] == 0:
        log_level = logging.INFO
    elif args['verbose'] == 1:
        log_level = logging.DEBUG
    else:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt="%y%m%d %H%M%S")

    results = convert_to_dict(args['pintron_genomic_file'],
                              args['pintron_align_file'],
                              args['pintron_introns_file'])

    with smart_open_out(args['output_json_file']) as fout:
        json.dump(results, fout)


if __name__ == "__main__":
    main()

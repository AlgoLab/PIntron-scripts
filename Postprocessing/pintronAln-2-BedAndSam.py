#!/usr/bin/env python

from __future__ import print_function

import argparse
import contextlib
import gzip
import json
import logging
import sys

import pybedtools


@contextlib.contextmanager
def smart_open_in(filename=None):
    if filename and filename != '-':
        fn = filename
        fo = open(filename, 'r')
        if filename.endswith(".gz"):
            fh = gzip.GzipFile(fn, 'rb', 9, fo)
        else:
            fh = fo
    else:
        fn = "<stdout>"
        fo = sys.stdin
        fh = sys.stdin

    try:
        yield fh
    finally:
        if fh is not fo:
            fh.close()
        if fo is not sys.stdin:
            fo.close()


def alignment_to_bed12(genomic_ref, alignment):
    blocks = alignment["blocks"]
    if genomic_ref["strand"] == "+":
        fstart = blocks[0]["genomic_absolute_start"]-1
        fend = blocks[-1]["genomic_absolute_end"]
        res = [genomic_ref["seqname"],
               fstart,
               fend,
               alignment["identifier"],
               "0",
               genomic_ref["strand"],
               fstart,
               fend,
               0,
               len(blocks),
               ",".join([str(block["genomic_absolute_end"] - block["genomic_absolute_start"]+1)
                         for block in blocks]),
               ",".join([str(block["genomic_absolute_start"] - fstart - 1)
                         for block in blocks])]
    else:
        fstart = blocks[-1]["genomic_absolute_end"] - 1
        fend = blocks[0]["genomic_absolute_start"]
        res = [genomic_ref["seqname"],
               fstart,
               fend,
               alignment["identifier"],
               "0",
               genomic_ref["strand"],
               fstart,
               fend,
               0,
               len(blocks),
               ",".join([str(block["genomic_absolute_start"] - block["genomic_absolute_end"] + 1)
                         for block in reversed(blocks)]),
               ",".join([str(block["genomic_absolute_end"] - fstart - 1)
                         for block in reversed(blocks)])]
    return res


def convert_pintron_align_to_bed12(genomic_block, alignments):
    logging.info("Starting conversion of %d alignments...", len(alignments))
    bed12 = []
    for alignment in alignments.values():
        bed12.append(alignment_to_bed12(genomic_block, alignment))
    logging.info("Conversion terminated.")
    return pybedtools.BedTool(bed12).sort()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--pintron-results-file',
                        default="-")
    parser.add_argument('-b', '--output-bed12-file',
                        type=argparse.FileType(mode='w'),
                        default=open("alignments.bed", "w"))
    parser.add_argument('-s', '--output-sam-file',
                        type=argparse.FileType(mode='w'),
                        default=open("alignments.sam", "w"))
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

    if args['output_bed12_file'] == args['output_sam_file']:
        sys.exit("The BED and SAM output files cannot be the same file!")

    logging.info("Reading PIntron results from file '%s'", args["pintron_results_file"])
    results = None
    with smart_open_in(args["pintron_results_file"]) as fjsonin:
        results = json.load(fjsonin)

    genomic_block = results["genome"]

    bed_desc = convert_pintron_align_to_bed12(genomic_block, results["alignments"])

    logging.info("Writing in BED format to file '%s'", args['output_bed12_file'].name)
    print(bed_desc, file=args['output_bed12_file'], end="")

    logging.info("Converting to SAM format...")
    sam_desc = bed_desc.to_bam(genome="hg38", bed12=True)
    logging.info("Writing in SAM format to file '%s'", args['output_sam_file'].name)
    print(sam_desc, file=args['output_sam_file'], end="")

    logging.info("Terminated.")


if __name__ == "__main__":
    main()

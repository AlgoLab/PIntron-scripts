#!/usr/bin/env python

from __future__ import print_function

import argparse
import logging
import sys

import pybedtools

class GenomicBlock:
    def __init__(self, genomic_header):
        logging.debug("Building block from string '%s'", genomic_header)
        parts = genomic_header.split(":")
        self.seqname = parts[0]
        self.start = int(parts[1])
        self.end = int(parts[2])
        self.strand = int(parts[3]+"1")
        self.strand = "+" if self.strand > 0 else "-"
        logging.debug("Built block '%s'", self)

    def __str__(self):
        return "{self.seqname}:{self.start}:{self.end}:{self.strand}".format(self=self)

    def to_bed(self):
        return (self.seqname, self.start-1, self.end, self.strand)

class SplicedAlignment:
    def __init__(self, identifier):
        logging.debug("New spliced alignment for sequence '%s'", identifier)
        self.identifier = identifier
        if self.identifier.startswith("/gb="):
            self.identifier = self.identifier[4:]
        self.blocks = []

    def __str__(self):
        return "{identifier} [{blocks}]".format(identifier=self.identifier,
                                                blocks=", ".join([str(block) for block in self.blocks]))

    def to_bed12(self, genomic_ref):
        assert self.blocks, "No blocks for the current alignment are available!"
        assert genomic_ref.strand == "+", "Calculations on reverse strand are not currently implemented!"
        fstart = genomic_ref.start + self.blocks[0].gstart-1
        fend = genomic_ref.start + self.blocks[-1].gend
        res = [genomic_ref.seqname,
               fstart,
               fend,
               self.identifier,
               "0",
               genomic_ref.strand,
               fstart,
               fend,
               0,
               len(self.blocks),
               ",".join([ str(block.gend-block.gstart+1) for block in self.blocks]),
               ",".join([ str(block.gstart-self.blocks[0].gstart) for block in self.blocks])]
        return res




class SplicedAlignmentBlock:
    def __init__(self, blockstr):
        parts = blockstr.split()
        (self.sstart, self.send,
         self.gstart, self.gend,
         self.sblock, self.gblock) = (int(parts[0]), int(parts[1]),
                                      int(parts[2]), int(parts[3]),
                                      parts[4], parts[5])

    def __str__(self):
        return "{self.gstart}-{self.gend}".format(self=self)

def pintron_alignments(alignments):
    curr_alignment = None
    for line in alignments:
        line = line.strip()
        if line.startswith(">"):      # New sequence/alignment
            if curr_alignment:
                yield curr_alignment
            curr_alignment = SplicedAlignment(line.lstrip(">").split()[0])
        elif line.startswith("#"):    # Annotation (currently ignored)
            pass
        elif line[0].isdigit():       # New block
            assert(curr_alignment)
            curr_alignment.blocks.append(SplicedAlignmentBlock(line))
    if curr_alignment:
        yield curr_alignment


def convert_pintron_align_to_bed12(infile, genomic_block):
    logging.info("Starting conversion of '%s'", infile.name)
    bed12 = []
    for alignment in pintron_alignments(infile):
        bed12.append(alignment.to_bed12(genomic_block))
    logging.info("Conversion terminated")
    return pybedtools.BedTool(bed12).sort()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pintron-align-file',
                        nargs='?',
                        type=argparse.FileType(mode='r'),
                        default=sys.stdin)
    parser.add_argument('-g', '--pintron-genomic-file',
                        nargs='?',
                        type=argparse.FileType(mode='r'),
                        default=open('genomic.txt'))
    parser.add_argument('-b', '--output-bed12-file',
                        nargs='?',
                        type=argparse.FileType(mode='w'),
                        default=sys.stdout)
    parser.add_argument('-s', '--output-sam-file',
                        nargs='?',
                        type=argparse.FileType(mode='w'),
                        default=sys.stdout)
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

    genomic_block = GenomicBlock(args['pintron_genomic_file'].readline().strip().lstrip(">"))

    bed_desc = convert_pintron_align_to_bed12(args['pintron-align-file'], genomic_block)

    logging.info("Writing in BED format to file '%s'", args['output_bed12_file'].name)
    print(bed_desc, file=args['output_bed12_file'])

    logging.info("Converting to SAM format...")
    sam_desc = bed_desc.to_bam(genome="hg38", bed12=True)
    logging.info("Writing in SAM format to file '%s'", args['output_sam_file'].name)
    print(sam_desc, file=args['output_sam_file'])




if __name__ == "__main__":
    main()

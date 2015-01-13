#!/usr/bin/python2.7

import pysam
from optparse import OptionParser
import os
import re
from progressbar import Bar, Timer, Percentage, ProgressBar

def main():
	parser = OptionParser(description = "Convert BAM file.",
                              usage = "%prog -B <aln_file.bam> -O <output_file.fa>")
        parser.add_option("-b",
                          metavar = "<aln_file.bam>",
                          help = "Alignment file in BAM format.")
        parser.add_option("-o",
                          metavar = "<output_file.fa>",
                          help = "Output file in FASTA format.")
        parser.add_option("-r",
                          metavar = "<chr:start-stop>",
                          help = "Region to be examined in the form: chr:start-stop. Ex. 1:100-200.",
                          default = "")
        (options, args) = parser.parse_args()
        in_bam_file = options.b
        out_fasta_file = options.o
        in_region = options.r

        if not (in_bam_file or out_file):
                print "Error: missing argument."
                return

        in_sam = pysam.AlignmentFile(in_bam_file, "rb")
        fetch_aln = in_sam.fetch()
        if in_region != "":
                regexp = re.compile(r'^(\w+):(\d+)-(\d+)')
                region = re.search(regexp, in_region)
                if region:
                        r_chr = region.group(1)
                        r_start = int(region.group(2))
                        r_stop = int(region.group(3))
                        fetch_aln = in_sam.fetch(r_chr, r_start, r_stop)
                        tot_fetch_aln = in_sam.count(r_chr, r_start, r_stop)
                        print "Parsing only read aligned to " + in_region
                else:
                        print "Error: wrong input region."
                        return
        else:
                print "Parsing all read alignments."
                tot_fetch_aln = in_sam.count()

        widgets = ['Processing: ', Percentage(), ' ', Bar(marker='=', left='[', right=']'), ' ', Timer()]
        bar = ProgressBar(widgets=widgets, maxval=tot_fetch_aln).start()

        with open(out_fasta_file, "w") as out_fasta:
                num_aln = 0
                for read in fetch_aln:
                        num_aln = num_aln + 1
                        bar.update(num_aln)
                        ref_name = in_sam.getrname(read.reference_id)
                        fasta_hdr = ">/id=" + read.query_name
                        if read.is_paired:
                                fasta_hdr += " /read=" + ("first" if read.is_read1 else "second")
                        fasta_hdr += " /reversed=" + ("yes" if read.is_reverse else "no")
                        fasta_hdr += " /ref_start=" + ref_name + ":" + str(read.reference_start)
                        fasta_hdr += " /ref_end=" + ref_name + ":" + str(read.reference_end)
                        #print fasta_hdr
                        out_fasta.write(fasta_hdr + "\n")
                        #print read.query_sequence
                        for b in read.get_blocks():
                                block_start = 0
                                block_stop = 0
                                for aln_pair in read.get_aligned_pairs():
                                        if aln_pair[1] == b[0]:
                                                block_start = aln_pair[0]
                                        if aln_pair[1] == b[1] - 1:
                                                block_stop = aln_pair[0]
                                # block_aln = "/query_block_start=" + str(block_start)
                                # block_aln += " /query_block_stop=" + str(block_stop)
                                # block_aln += " /ref_block_start=" + ref_name + ":" + str(b[0])
                                # block_aln += " /ref_block_stop=" + ref_name + ":" + str(b[1])
                                block_aln = str(block_start) + " " + str(block_stop)
                                block_aln += " " + ref_name + ":" + str(b[0])
                                block_aln += " " + ref_name + ":" + str(b[1])
                                #print block_aln
                                out_fasta.write(block_aln + "\n")
                                #print read.query_sequence[block_start:block_stop + 1]
                                out_fasta.write(read.query_sequence[block_start:block_stop + 1] + "\n")
                                #print read.get_aligned_pairs()
                out_fasta.close()
                bar.finish()
                print "Num. Processed Alignments: {0}".format(num_aln)
        in_sam.close()

if __name__ == '__main__':
        main()

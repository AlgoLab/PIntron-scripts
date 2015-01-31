#!/usr/bin/python2.7

import pysam
from optparse import OptionParser
import os
import re
from progressbar import Bar, Timer, Percentage, ProgressBar

def main():
        parser = OptionParser(description = "Split aligned read from BAM file "\
                              "according to the annotation file.",
                              usage = "%prog -b <aln_file.bam> -a <annotation_file.csv> "\
                              "[-r <chr:start-stop> | -g <gene_name>]")
        parser.add_option("-b",
                          metavar = "<aln_file.bam>",
                          help = "Alignment file in BAM format.")
        parser.add_option("-r",
                          metavar = "<chr:start-stop>",
                          help = "Region to be examined in the form: chr:start-stop. Ex. 1:100-200.",
                          default = "")
        parser.add_option("-g",
                          metavar = "<gene_name>",
                          help = "Gene name.",
                          default = "")
	parser.add_option("-a",
                          metavar = "<annotation_file.csv>",
                          help = "Gene annotation file in CSV format.")
        (options, args) = parser.parse_args()
        in_bam_file = options.b
        in_region = options.r
        in_gene = options.g
	in_annot_file = options.a

        if not (in_bam_file and in_annot_file):
		print "Error: missing argument."
		return

        if (in_region != "" and in_gene != ""):
                print "Error: specify only one option between region (-r) and gene name (-g)."
                return

        regexp_reg = re.compile(r'^(\w+):(\d+)-(\d+)')
	regexp_annot = re.compile(r'^([\dXY]+)\t(\d+)\t(\d+)\t([+\-])' \
		'\tgene_id=\"(\w+)\"\tgene_version=\"(\d+)\"\tgene_name=\"([\w\-\.]+)\"')

	in_annot = open(in_annot_file, "r")

        in_sam = pysam.AlignmentFile(in_bam_file, "rb")
        fetch_aln = in_sam.fetch()
	tot_fetch_aln = 0
	count = 0
	match_elem = {}

        if(in_region != ""):
                region = re.search(regexp_reg, in_region)
                if(region):
                        in_chr = region.group(1)
                        in_start = int(region.group(2))
                        in_stop = int(region.group(3))
                else:
                        print "Error: wrong input region."
                        return

	for a in in_annot:
		ga = re.search(regexp_annot, a)
		if(ga):
			count = count + 1
			chr = ga.group(1)
			start = int(ga.group(2))
			stop = int(ga.group(3))
			strand = ga.group(4)
			gene_id = str(ga.group(5))
			gene_version = int(ga.group(6))
			gene_name = str(ga.group(7))
			insert = False
			el = {'chr' : chr,
				'start' : start,
				'stop' : stop,
				'strand' : strand,
				'gene_id' : gene_id,
				'gene_version' : gene_version,
				'gene_name' : gene_name}
			if(in_gene != ""):
				if(in_gene == gene_name or in_gene == gene_id):
					insert = True
			elif(in_region != ""):
				if(in_chr == chr and in_start <= start and in_stop >= stop):
					insert = True
			else:
				insert = True
			if(insert):
				if not (match_elem.has_key(gene_name)):
					match_elem[gene_name] = el
				elif (match_elem[gene_name]['gene_version'] < gene_version):
					match_elem[gene_name] = el
				else:
					print "WARNING: Ducplicated gene name";
				#print a
		else:
			print "Error in parsing annotation file."
			print a
			return
	print "Parsed {0} annotated genes.".format(count)
        print "Num. of retrieved genes: {0}.".format(len(match_elem))
        for k in match_elem.keys():
                print "Creating gene {0}.".format(k)
                r_chr = match_elem[k]['chr']
                r_start = match_elem[k]['start'] - 1000
                r_stop = match_elem[k]['stop'] + 1000
                fetch_aln = in_sam.fetch(r_chr, r_start, r_stop)
                tot_fetch_aln = in_sam.count(r_chr, r_start, r_stop)
                if(tot_fetch_aln == 0):
                        print "No valid alignments found."
                        continue
                print "Num. aligned reads: {0}".format(tot_fetch_aln)

                widgets = ['Processing: ', Percentage(), 
                           ' ', Bar(marker='=', left='[', right=']'),
                           ' ', Timer()]
                bar = ProgressBar(widgets=widgets, maxval=tot_fetch_aln).start()
                if not (os.path.exists(k)):
                        os.mkdir(k)
                with open(k + "/" + k + ".fa", "w") as out_fasta:
                        num_aln = 0
                        for read in fetch_aln:
                                num_aln = num_aln + 1
                                bar.update(num_aln)
                                ref_name = in_sam.getrname(read.reference_id)
                                fasta_hdr = ">/gb=" + read.query_name
                                if read.is_paired:
                                        fasta_hdr += ("_R1" if read.is_read1 else "_R2")
                                fasta_hdr += " /clone_end=3'" + " /reversed="
                                fasta_hdr += ("yes" if read.is_reverse else "no")
                                fasta_hdr += " /ref_start=" + ref_name
                                fasta_hdr += ":" + str(read.reference_start)
                                fasta_hdr += " /ref_end=" + ref_name
                                fasta_hdr += ":" + str(read.reference_end)
                                #print fasta_hdr
                                out_fasta.write(fasta_hdr + "\n")
                                #print read.query_sequence
                                out_fasta.write(read.query_sequence + "\n")
                        out_fasta.close()
                        bar.finish()
                        print "Num. Processed Alignments: {0}".format(num_aln)
        in_sam.close()

if __name__ == '__main__':
        main()

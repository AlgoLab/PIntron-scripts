#!/usr/bin/python2.7

import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
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
        parser.add_option("-f",
                          metavar = "<fasta_dir>",
                          help = "Directory containing FASTA (.gz) files of the chromosomes.")
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
        in_fasta_dir = options.f
        in_fasta_prefix_name = "Homo_sapiens.GRCh38.dna.chromosome."
        in_gene = options.g
	in_annot_file = options.a

        if not (in_bam_file and in_annot_file and in_fasta_dir):
		print "Error: missing argument."
		return

        if (in_region != "" and in_gene != ""):
                print "Error: specify only one option between region (-r) and gene name (-g)."
                return

        regexp_reg = re.compile(r'^(\w+):(\d+)-(\d+)')
	regexp_annot = re.compile(r'^([\dXY]+)\t(\d+)\t(\d+)\t([+\-])' \
		'\tgene_id=\"(\w+)\"\tgene_version=\"(\d+)\"\tgene_name=\"([\w\-\.]+)\"')

        if not (os.path.exists(in_fasta_dir)):
                print "Error: directory {0} not found.".format(in_fasta_dir)

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
			strand = str(ga.group(4))
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
                print ""
                print "Creating gene {0}.".format(k)
                r_chr = match_elem[k]['chr']
                r_start = match_elem[k]['start'] - 1000
                r_stop = match_elem[k]['stop'] + 1000
                r_strand = match_elem[k]['strand']
                fetch_aln = in_sam.fetch(r_chr, r_start, r_stop)
                tot_fetch_aln = in_sam.count(r_chr, r_start, r_stop)
                if(tot_fetch_aln == 0):
                        print "No valid alignments found."
                        continue

                if not (os.path.exists(k)):
                        os.mkdir(k)

                #Compute reads
                widgets = ['Processing: ', Percentage(),
                           ' ', Bar(marker='=', left='[', right=']'),
                           ' ', Timer()]
                bar = ProgressBar(widgets=widgets, maxval=tot_fetch_aln).start()
                with open(k + "/" + k + ".fa", "w") as out_fasta:
                        num_proc_seq = 0
                        num_valid_seq = 0
                        for read in fetch_aln:
                                num_proc_seq = num_proc_seq + 1
                                bar.update(num_proc_seq)
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
                                #print read.query_sequence
                                if not (read.is_paired):
                                        num_valid_seq = num_valid_seq + 1
                                        out_fasta.write(fasta_hdr + "\n")
                                        out_fasta.write(read.query_sequence + "\n")
                                elif not (read.mate_is_unmapped):
                                        num_valid_seq = num_valid_seq + 1
                                        out_fasta.write(fasta_hdr + "\n")
                                        out_fasta.write(read.query_sequence + "\n")
                        out_fasta.close()
                        bar.finish()
                        print "Num. Processed Sequences: {0}".format(num_proc_seq)
                        print "Num. Valid Sequences: {0}".format(num_valid_seq)

                #Compute genomics
                print "Cutting genomic sequence."
                found = True
                results = []
                seq_name = ""
                if not(r_chr[:3] == "chr"):
                        seq_name += "chr"
                fasta_seq_name = in_fasta_dir + "/" + in_fasta_prefix_name + r_chr + ".fa.gz"
                if not os.path.isfile(fasta_seq_name):
                        print "Error: file {0} does not exists.".format(fasta_seq_name)
                        return
                seq_handle = gzip.open(fasta_seq_name, "r")
                for sequence in SeqIO.parse(seq_handle, "fasta"):
                        if (sequence.id == r_chr):
                                found = True
                                sub_s = Seq(str(sequence.seq[r_start:r_stop]), generic_dna)
                                seq_id = "{0}{1}:{2}:{3}:".format(seq_name,r_chr,r_start,r_stop)
                                if(r_strand == "-"):
                                        sub_s = sub_s.reverse_complement()
                                        seq_id += "-1"
                                        #descr += " ReverseComplemented"
                                else:
                                        seq_id += "+1"
                                seqrec = SeqRecord(sub_s)
                                seqrec.id = str(seq_id)
                                #descr += " Length={0}bp.".format(len(sub_s))
                                seqrec.description = ""
                                results.append(seqrec)
                                #print seqrec.seq
                                print "Cut sequence of {0}bp.".format(len(sub_s))
                if not (found):
                        print "No sequence {0} found in {1}.".format(chr, in_fasta_file)
                else:
                        out_genomic = open(k + "/" + "genomic.txt", "w")
                        SeqIO.write(results, out_genomic, "fasta")
                        out_genomic.close()

        in_sam.close()

if __name__ == '__main__':
        main()

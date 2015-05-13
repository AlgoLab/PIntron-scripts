#!/usr/bin/python2.7

import logging
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam
from optparse import OptionParser
import os
import re
import sys
from progressbar import Bar, Timer, Percentage, ProgressBar

def main():
        parser = OptionParser(description = "Split aligned read from BAM file "\
                              "according to the annotation file.",
                              usage = "%prog -b <aln_file.bam> -a <annotation_file.csv> "\
                              "-f <fasta_dir> [-r <chr:start-stop> | -g <gene_name>]")
        parser.add_option("-b",
                          metavar = "<aln_file.bam>",
                          help = "Alignment file in BAM format.")
	parser.add_option("-a",
                          metavar = "<annotation_file.csv>",
                          help = "Gene annotation file in CSV format.")
        parser.add_option("-f",
                          metavar = "<fasta_dir>",
                          help = "Directory containing FASTA (.gz) files of the chromosomes.")
        parser.add_option("-r",
                          metavar = "<chr:start-stop>",
                          help = "Region to be examined in the form: chr:start-stop. Ex. 1:100-200.",
                          default = "")
        parser.add_option("-g",
                          metavar = "<gene_name>",
                          help = "Gene name.",
                          default = "")
        parser.add_option("-o",
                          metavar = "<output-dir>",
                          help = "Output (root) directory.",
                          default = ".")
        parser.add_option('-v',
                        help='increase output verbosity',
                        action='count', default=0)
        (options, args) = parser.parse_args()
        in_bam_file = options.b
        in_region = options.r
        in_fasta_dir = options.f
        in_fasta_prefix_name = "Homo_sapiens.GRCh38.dna.chromosome."
        in_gene = options.g
	in_annot_file = options.a
        out_root_dir = options.o

        if options.v == 0:
                log_level = logging.INFO
        elif options.v == 1:
                log_level = logging.DEBUG
        else:
                log_level = logging.DEBUG

        logging.basicConfig(level=log_level,
                            format='%(levelname)-8s [%(asctime)s]  %(message)s',
                            datefmt="%y%m%d %H%M%S")

        if not (in_bam_file and in_annot_file and in_fasta_dir):
		logging.error("Missing input argument(s).")
		sys.exit(1)

        if (in_region != "" and in_gene != ""):
                logging.error("Specify only one option between region (-r) and gene name (-g).")
                sys.exit(1)

        if not (os.path.exists(out_root_dir)):
                logging.error("Output dir not found.")
                sys.exit(1)

        regexp_reg = re.compile(r'^(\w+):(\d+)-(\d+)')
	regexp_annot = re.compile(r'^([\dXY]+)\t(\d+)\t(\d+)\t([+\-])' \
		'\tgene_id=\"(\w+)\"\tgene_version=\"(\d+)\"\tgene_name=\"([\w\-\.]+)\"')

        if not (os.path.exists(in_fasta_dir)):
                logging.error("Directory " + in_fasta_dir + " not found.")
                sys.exit(1)

        logging.info("splitBAM: Program Started")

	in_annot = open(in_annot_file, "r")

        in_sam = pysam.AlignmentFile(in_bam_file, "rb")
        fetch_aln = in_sam.fetch()
	tot_fetch_aln = 0
	count = 0
	match_elem = {}

        if(in_region != ""):
                region = re.search(regexp_reg, in_region)
                if (region):
                        in_chr = region.group(1)
                        in_start = int(region.group(2))
                        in_stop = int(region.group(3))
                else:
                        logging.error("Wrong input region.")
                        sys.exit(1)
                if (in_start > in_stop):
                        logging.error("Wrong input region.")
                        sys.exit(1)

	for a in in_annot:
		ga = re.search(regexp_annot, a)
		if (ga):
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
					logging.warn("Ducplicated gene name")
				logging.debug(a)
		else:
			logging.error("Error in parsing annotation file.")
			logging.error(a)
			sys.exit(1)
        logging.info("Parsed " + str(count) + " annotated genes.")
        logging.info("Num. of retrieved genes: " + str(len(match_elem)))
        for k in match_elem.keys():
                logging.info("")
                logging.info("Creating gene " + k)
                r_chr = match_elem[k]['chr']
                r_start = match_elem[k]['start'] - 1000
                r_stop = match_elem[k]['stop'] + 1000
                r_strand = match_elem[k]['strand']
                fetch_aln = in_sam.fetch(r_chr, r_start, r_stop)
                tot_fetch_aln = in_sam.count(r_chr, r_start, r_stop)
                if(tot_fetch_aln == 0 or tot_fetch_aln < 10):
                        logging.warn("No valid alignments found.")
                        continue

                out_dir = out_root_dir + "/" + k
                if not (os.path.exists(out_dir)):
                        os.mkdir(out_dir)

                #Compute reads
                widgets = ['Processing: ', Percentage(),
                           ' ', Bar(marker='=', left='[', right=']'),
                           ' ', Timer()]
                bar = ProgressBar(widgets=widgets, maxval=tot_fetch_aln).start()
                num_proc_seq = 0
                num_valid_seq = 0
                num_disc_seq = 0
                valid_id = set()
                valid = []
                discarded = []
                for read in fetch_aln:
                        num_proc_seq = num_proc_seq + 1
                        bar.update(num_proc_seq)
                        ref_name = in_sam.getrname(read.reference_id)
                        read_name = read.query_name
                        if read.is_paired:
                                read_name += ("_R1" if read.is_read1 else "_R2")
                        fasta_hdr = "/gb=" + read_name
                        fasta_hdr += " /clone_end=3'" + " /reversed="
                        fasta_hdr += ("yes" if read.is_reverse else "no")
                        fasta_hdr += " /ref_start=" + ref_name
                        fasta_hdr += ":" + str(read.reference_start)
                        fasta_hdr += " /ref_end=" + ref_name
                        fasta_hdr += ":" + str(read.reference_end)
                        record = SeqRecord(Seq(read.query_sequence, generic_dna))
                        record.id = fasta_hdr
                        record.description = ""
                        logging.debug(fasta_hdr)
                        logging.debug(read.query_sequence)
                        if not (read.is_paired):
                                if read_name not in valid_id:
                                        num_valid_seq = num_valid_seq + 1
                                        valid.append(record)
                                        valid_id.add(read_name)
                        else:
                                if not (read.mate_is_unmapped):
                                        if (read.reference_id == read.next_reference_id and
                                            read.next_reference_start >= r_start and
                                            read.next_reference_start <= r_stop):
                                                if read_name not in valid_id:
                                                        num_valid_seq = num_valid_seq + 1
                                                        valid.append(record)
                                                        valid_id.add(read_name)
                                        else:
                                                num_disc_seq = num_disc_seq + 1
                                                discarded.append(record)
                                else:
                                        num_disc_seq = num_disc_seq + 1
                                        discarded.append(record)
                bar.finish()
                out_fasta = open(out_dir + "/" + k + ".fa", "w")
                SeqIO.write(valid, out_fasta, "fasta")
                out_fasta.close()
                out_dis = open(out_dir + "/" + "discarded.fa", "w")
                SeqIO.write(discarded, out_dis, "fasta")
                out_dis.close()
                logging.info("Num. Processed Sequences: " + str(num_proc_seq))
                logging.info("Num. Valid Sequences: " + str(num_valid_seq))
                logging.info("Num. Discarded Sequences: " + str(num_disc_seq))

                #Compute genomics
                logging.info("Cutting genomic sequence.")
                found = True
                sequences = []
                seq_name = ""
                if not(r_chr[:3] == "chr"):
                        seq_name += "chr"
                fasta_seq_name = in_fasta_dir + "/" + in_fasta_prefix_name + r_chr + ".fa.gz"
                if not os.path.isfile(fasta_seq_name):
                        logging.error("File " + fasta_seq_name + " not found.")
                        sys.exit(1)
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
                                sequences.append(seqrec)
                                logging.debug(seqrec.seq)
                                logging.info("Cut sequence of " + str(len(sub_s)) + "bp.")
                if not (found):
                        logging.warn("No sequence " + chr + " found in " + in_fasta_file)
                else:
                        out_genomic = open(out_dir + "/" + "genomic.txt", "w")
                        SeqIO.write(sequences, out_genomic, "fasta")
                        out_genomic.close()

        in_sam.close()
        logging.info("splitBAM: Program Completed")

if __name__ == '__main__':
        main()

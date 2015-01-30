#!/usr/bin/python2.7

import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from optparse import OptionParser
import re

def main():
        parser = OptionParser(description = "Cut geomic sequence.",
                              usage = "%prog -F <fasta_file.fa[.gz]> -O <output_file.fa> -R <chr:start-stop>")
        parser.add_option("-F",
                          metavar = "<fasta_file.fa[.gz]>",
                          help = "Fasta original file containing the full sequence. Optionally in GZIP format.")
        parser.add_option("-O",
                          metavar = "<output_file.fa>",
                          help = "Output file in FASTA format.")
        parser.add_option("-R",
                          metavar = "<chr:start-stop>",
                          help = "Region to be cut in the form: chr:start-stop. Ex. 1:100-200.")
	parser.add_option("-D",
                          metavar = "<sequence_description>",
                          help = "Description of the sequence added in the FASTA header.",
			  default="")
	parser.add_option("--rev_comp",
                          action="store_true",
			  dest="flag",
                          help = "Description of the sequence added in the FASTA header.",
                          default=False)
        (options, args) = parser.parse_args()
        in_fasta_file = options.F
        out_fasta_file = options.O
        cut_region = options.R
	descr = options.D
	revcomp = options.flag

	regexp = re.compile(r'^(\w+):(\d+)-(\d+)')
	region = re.search(regexp, cut_region)
	if region:
		chr = region.group(1)
		start = int(region.group(2))
		stop = int(region.group(3))
		print "Cut region of sequence {0} from {1} to {2}.".format(chr, start, stop)
	else:
		print "Error: wrong input region."
		return

	handle = gzip.open(in_fasta_file, "r")
	found = False
	results = []
	for sequence in SeqIO.parse(handle, "fasta"):
		if (sequence.id == chr):
			found = True
			sub_s = Seq(str(sequence.seq[start:stop]), generic_dna)
			seq_id = "chr{0}:{1}:{2}:".format(chr,start,stop)
			if(revcomp):
				sub_s = sub_s.reverse_complement()
				seq_id += "-1"
				#descr += " ReverseComplemented"
			else:
				seq_id += "+1"
			seqrec = SeqRecord(sub_s)
			seqrec.id = str(seq_id)
			#descr += " Length={0}bp.".format(len(sub_s))
			seqrec.description = str(descr)
			results.append(seqrec)
			#print seqrec.seq
			print "Cut sequence of {0}bp.".format(len(sub_s))

	if not (found):
		print "No sequence {0} found in {1}.".format(chr, in_fasta_file)
	else:
		out_fasta = open(out_fasta_file, "w")
		SeqIO.write(results, out_fasta, "fasta")
		out_fasta.close()


if __name__ == '__main__':                                                                            
        main() 

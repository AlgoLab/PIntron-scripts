## Goals

To process the files that have been output by PIntron, so that the result can be visualized by pintronweb.

There are two steps:
1. to process the introns and the alignments, producing a list of all exons detected. *cutfiller* is the program for this phase.
2. to process the exons, producing the gene structure, *cutassemble* is the program for this phase.

## cutfiller

This script reconstructs gene putative exons from alignments of RNA-seq data on a genome.

### Input data

- file ```out-after-intron-agree.txt``` of the alignments produced by PIntron (option -a or --alignment-file; default: none)
- file ```predicted-introns.txt``` of the introns produced by PIntron (option -p or --pintron-output; default: none)
- maximum accepted intron length. Introns exceeding this value will be discarded (option -l or --max-intron-length; default: 15000)

As an alternative to the above files, cutfiller can take in input a SAM file containing the RNA-seq alignments (option -s or --samfile; default: none)

Each alignment in the input file out-after-intron-agree.txt is represented in the following format:
	
	><header>
	#polya=<polya_flag>
	#polyad=<polyadenilation_flag>
	<aligned_factor>+
	
	<aligned_factor> is a record of the following six space-separated fields:
		- start of the aligned read factor
		- end of the aligned read factor
		- start of the aligned genome factor (relative to the input genome of PIntron)
		- end of the aligned genome factor (relative to the input genome of PIntron)
		- nt sequence of the read factor
		- nt sequence of the genome factor
		
	In <header> the GenBank ID of the aligned read is specified after the string "/gb="
	
Each intron the input file predicted-introns.txt is represented as a record of the following twenty tab-separated fields:

		- start on the genome (input of PIntron)
		- end on the genome (input of PIntron)
		- start on the chromosome
		- end on the chromosome
		- intron length
		- number of reads supporting the intron
		- comma-separated list of the reads (GenBank IDs) supporting the intron
		- mean alignment error at the donor (5') site
		- mean alignment error at the acceptor (3') site
		- donor score
		- acceptor score
		- BPS score
		- BPS position
		- intron type (U2, U12, or unclassified)
		- intron pattern
		– repeat sequence
		– donor exon suffix
		– intron prefix
		– intron suffix
		– acceptor exon prefix


###Output data

A list of exons, encoded by:
		the exon start on the genome (input of PIntron)
		the exon end on the genome (input of PIntron)
		a flag indicating that the exon is complete (flag=0), is uncomplete on the left (flag=1) or is uncomplete on the right
		(flag=2)

The output format is still undecided, since it is mainly for internal uses to interact with ```cutassemble```.



Log information is produced in standard error.
	
###Options
* input files
* output file (default: stdout)
* log file (default: stderr)

###Cutassemble

Analyzes the results of ```cutfiller``` and produces a json file that is suitable for the visualization part. It should conform with the format of PIntron output file.

###Options
* input files
* output file (default: stdout)
* log file (default: stderr)

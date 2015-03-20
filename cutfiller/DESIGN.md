#  Algorithm of cutfiller #

Input:
- the set I of the introns on the genome
- the alignments of the RNA-seq reads on the genome
		
Output:
- the set of the putative exons

## Definitions ##

A splice site s=(l,l+1) can be one of the following things:
- an exon-intron boundary (l is the end of the exon at 5' and l+1 is the start of the intron at 3')
- an intron-exon boundary (l is the end of the intron at 5' and l+1 is the start of the exon at 3')
	
A splice site might also be both an exon-intron and an intron-exon boundary.
	
Types of ss:
- "5-type" if the ss is an exon-intron boundary
- "3-type" if the ss is an intron-exon boundary
		
An intron (s,e), where s and e are the start and the end respectively, gives the two following ss:
- a "5-type" splice site (s-1,s)
- a "3-type" splice site (e,e+1)

A ss (l,l+1) that is both "5-type" and "3-type" is given by two different introns (s,l) and (l+1,e).

## Steps of the algorithm ##

- STEP1: compute from the intron set I the set S5 of the "5-type" splice sites and the set S3 of the "3-type" splice sites.
Notice that S5 and S3 are not necessarily disjoint. Each splice site in the two sets are labeled as "unused".
	
- STEP2: for each pair (s1, s2) in the set (S3 X S5) such that s2=(l2,l2+1) is after s1=(l1,l1+1) on the genome, then the
mapping coverage from l1+1 to l2 is checked. If this region is completely covered by the read mapping,
then (l1+1, l2) is a complete putative exon. In that case, s1 and s2 are labeled as "used".
	
- STEP3: for each "unused" splice site (l,l+1) in S5, an uncomplete (on the left) exon (s,l) is produced, if the region
from s to l is covered by the the read mapping (maximal region).
	
- STEP4: for each "unused" splice site (l,l+1) in S3, an uncomplete (on the right) exon (l+1,e) is produced, if the region
from l+1 to e is covered by the the read mapping (maximal region).

- STEP5: for each "used" splice site (l,l+1) in S5, an uncomplete (on the left) exon (s,l) is produced, if (i) the region
from s to l is covered by the the read mapping (maximal region), (ii) s is less than the minimum start sm of the
complete exons ending in l (and computed by STEP1), and (iii) (sm-s) is over a given threshold.
		
- STEP6: for each "used" splice site (l,l+1) in S3, an uncomplete (on the right) exon (l+1,e) is produced, if (i) the region
from l+1 to e is covered by the the read mapping (maximal region), (ii) e is greater than the maximum end em of the
complete exons starting in l+1 (and computed by STEP1), and (iii) (e-em) is over a given threshold.

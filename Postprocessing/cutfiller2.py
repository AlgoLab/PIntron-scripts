#!/usr/bin/env python2.7

import logging
import sys
import os
import argparse
import numpy
import pysam
import json

def is_unique(region, offset):
    """
    Analyze region and return true if it represents a single,
    contiguous region.
    """
    if len(region) == 0:
        return -1
    else:
        uniq = True
        i = 0
        s = ""
        while (region[i] == 0):
            i += 1
        s = str(i + offset)
        is_zero = False
        begin = i
        while (i < len(region)):
            if(region[ i ] == 0):
                if(not is_zero):
                    s = s + "--" + str(i + offset - 1) + "//"
                    is_zero = True
                uniq = False
            else:
                if(is_zero):
                    s = s + str(i + offset)
                    is_zero = False
            i += 1
        i = -1
        while(region[ i ] == 0):
            i -= 1
        s = s + "--" + str(len(region) + i + offset)
        logging.debug("==> COVERAGE: {0}".format(s))
        end = len(region) + i
        if(not uniq):
            return 0
        else:
            reg_min = region[ begin : end ].min()
            reg_max = region[ begin : end ].max()
            logging.debug("Min Cov: {0} -- Max Cov: {1}".format(reg_min, reg_max))
            if(reg_min < reg_max / 10):
                return 0
            else:
                mean_cov = region[ begin : end ].sum() / (end - begin)
                logging.debug("Mean cov: {0}".format(int(round(mean_cov))))
                return int(round(mean_cov))


def grow_right (region, offset, max_gap):
    """
    Compute the biggest contiguous region that starts from the first
    element of the region (if any).
    """
    if len(region) == 0:
        return []
    else:
        exon = [ 0, 0, 0 ]
        i = 0
        begin = i
        max_b = region[ 0 ]
        while(i < max_gap):
            if(region[ i ] > max_b):
                max_b = region[ i ]
                begin = i
            i += 1
        exon[ 0 ] = begin + offset
        stop = False
        while(i < len(region) and region[ i ] != 0 and not stop):
            reg_min = region[ begin : i ].min()
            reg_max = region[ begin : i ].max()
            if(reg_min >= reg_max / 10):
                end = i
            else:
                stop = True
            i += 1
        exon[ 1 ] = i + offset
        if(i == len(region)):
            exon[ 1 ] -= 1
        
        logging.debug("Min Cov: {0} -- Max Cov: {1}".format(region[ begin : end ].min(),
                                                            region[ begin : end ].max()))
        mean_cov = region[ begin : end ].sum() / (end - begin)
        logging.debug("Mean cov: {0}".format(int(round(mean_cov))))
        exon[ 2 ] = int(round(mean_cov))
        numpy.set_printoptions(threshold='nan')
        logging.debug(region[ begin : end ])
        return exon

def grow_left(region, offset, max_gap):
    """
    Compute the biggest contiguous region that ends from the last
    element of the region (if any).
    """
    if len(region) == 0:
        return []
    else:
        exon = [ 0, 0, 0 ]
        i = 1
        end = len(region) - 1
        max_e = region[ end ]
        while(i < max_gap):
            if(region[ len (region) - i ] > max_e):
                end = len(region) - i
                max_e = region[ end ]
            i += 1
        exon[ 1 ] = end + offset + 1
        stop = False
        while(i < len(region) and region[ i ] != 0 and not stop):
            reg_min = region[ len(region) - i : end ].min()
            reg_max = region[ len(region) - i : end ].max()
            if(reg_min >= reg_max / 10):
                begin = len(region) - i
            else:
                stop = True
            i += 1
        exon[ 0 ] = begin + offset
        logging.debug("Min Cov: {0} -- Max Cov: {1}".format(region[ begin : end ].min(),
                                                            region[ begin : end ].max()))
        mean_cov = region[ begin : end ].sum() / (end - begin)
        logging.debug("Mean cov: {0}".format(int(round(mean_cov))))
        exon[ 2 ] = int(round(mean_cov))
        numpy.set_printoptions(threshold='nan')
        logging.debug(region[ begin : end ])
        return exon


def main():
    parser = argparse.ArgumentParser(prog = "cutfiller",
                                      description = "Rebuild exons from PIntron output.",
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-j', '--json-output', help = "Full output in JSON format.",
                        required = True, dest = 'json_file')
    parser.add_argument('-l', '--max-intron-length',
                         help = 'Max intron length accepted. Introns that exceed this value will be discarded',
                         required = False, dest = 'maxIntron', type = int, default = 100000)
    parser.add_argument('-g', '--gap-ends', help = 'Gap threshold at exon ends for their computation.',
                         required = False, dest = 'gapEnds', type = int, default = 3)
    parser.add_argument('-v', '--verbose',
                        help='increase output verbosity',
                        action='count', default=0)
    args = parser.parse_args()

    if args.verbose == 0:
        log_level = logging.INFO
    elif args.verbose == 1:
        log_level = logging.DEBUG
    else:
        log_level = logging.DEBUG
        
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',
                        datefmt="%y%m%d %H%M%S")
    
    if not args.json_file:
        logging.error('No JSON file given.\nAborting...')
        sys.exit(1)

    logging.info("==> OPENING JSON FILE")
    jdata = {}
    if args.json_file.endswith(".gz"):
       with gzip.open(args.json_file, 'rb') as jfile:
           jdata = json.loads(jfile.read().decode("ascii"))
    else:
        with open(args.json_file, 'r') as jfile:
            jdata = json.load(jfile)

    logging.info("==> PROCESSING GENOME DATA")
    genome = jdata['genome']
    seq_name = genome['seqname']
    seq_start = genome['start']
    seq_end = genome['end']
    seq_strand = genome['strand']
    logging.debug("Name: {0}, Start: {1}, End: {2}, Strand: {3}".format(seq_name, seq_start, seq_end, seq_strand))

    logging.info("==> PROCESSING ALIGNMENT DATA")
    aln = jdata['alignments']
    align_index = 0
    for a in aln.keys():
        logging.debug("Aln Id: " + aln[a]['identifier'])
        for b in aln[a]['blocks']:
            logging.debug(str(b['genomic_relative_start']) + " " + str(b['genomic_relative_end']))
            logging.debug(str(b['genomic_absolute_start']) + " " + str(b['genomic_absolute_end']))
            r_start = b['genomic_relative_start']
            r_end = b['genomic_relative_end']
            abs_start = b['genomic_absolute_start']
            abs_end = b['genomic_absolute_end']
            if(align_index == 0):
                min_aln_pos = r_start
                max_aln_pos = r_end
            if r_start < min_aln_pos:
                min_aln_pos = r_start
            if r_end > max_aln_pos:
                max_aln_pos = r_end
            align_index += 1
    cov_reads = numpy.zeros(max_aln_pos - min_aln_pos + 1, dtype=numpy.uint32)
    align_index = 0
    for a in aln.keys():
        for b in aln[a]['blocks']:
            r_start = b['genomic_relative_start']
            r_end = b['genomic_relative_end']
            abs_start = b['genomic_absolute_start']
            abs_end = b['genomic_absolute_end']
            cov_reads[ r_start - min_aln_pos : r_end - min_aln_pos + 1 ] += 1
            align_index += 1
    logging.info("==> PROCESSED " + str(align_index) + " ALIGNMENTS")
    logging.debug("Min: {0} -- Max: {1}".format(min_aln_pos, max_aln_pos))
    logging.debug("Size: " + str(len(cov_reads)))

    five_p_sites = []
    five_p_sites_no_exons = []
    three_p_sites = []
    three_p_sites_no_exons = []

    logging.info("==> PROCESSING INTRON DATA")
    int_num = 1
    introns = jdata['introns']
    for i in introns:
        region_end = i['relative_start'] - 1
        region_begin = i['relative_end'] + 1
        abs_start = i['absolute_start']
        abs_end = i['absolute_end']
        int_num += 1
        if abs(region_begin - region_end) <  args.maxIntron:
            if region_end not in five_p_sites:
                five_p_sites.append(region_end)
                five_p_sites_no_exons.append (region_end)
            if region_begin not in three_p_sites:
                three_p_sites.append(region_begin)
                three_p_sites_no_exons.append(region_begin)
        else:
            logging.info("==> INTRON TOO LARGE: [{0}, {1}]".format(region_begin, region_end))

    three_p_sites = sorted(list(set(three_p_sites)))
    five_p_sites = sorted(list(set(five_p_sites)))
    five_p_sites_no_exons = sorted(list(set(five_p_sites_no_exons)))
    three_p_sites_no_exons = sorted(list(set(three_p_sites_no_exons)))

    logging.info("==> PROCESSING 3' SPLICE SITES")
    exon_id = 1
    # Analyze complete regions for every splice site
    for site1 in three_p_sites:
        logging.info("==> COMPUTING NUMBER OF COMPLETE REGIONS FOR STARTING SITE {0}".format(site1))
        complete_regions_for_site1 = []
        should_continue = True
        for site2 in five_p_sites:
            if should_continue and site1 < site2:
                reg = cov_reads[ site1 - min_aln_pos : site2 - min_aln_pos ]
                if(reg.sum() > 0 and
                    cov_reads[ site1 - min_aln_pos - args.gapEnds : site1 - min_aln_pos + args.gapEnds ].sum() > 0 and
                    cov_reads[ site2 - min_aln_pos - args.gapEnds : site2 - min_aln_pos + args.gapEnds ].sum() > 0
               ):
                    logging.info("==> INTERVAL: {0}--{1}".format(site1, site2))
                    cov = is_unique(reg, site1)
                    if (cov > 0):
                        if site1 in three_p_sites_no_exons:
                            three_p_sites_no_exons = filter(lambda x: x != site1, three_p_sites_no_exons)
                        if site2 in five_p_sites_no_exons:
                            five_p_sites_no_exons = filter(lambda x: x != site2, five_p_sites_no_exons)
                        complete_regions_for_site1.append((site1, site2, cov))
                    else:
                        should_continue = False
                else:
                    logging.info("==> INTERVAL: {0}--{1}".format(site1, site2))
                    if(cov_reads[ site1 - min_aln_pos : site2 - min_aln_pos ].sum() > 0):
                        logging.info("==> EXCEEDING BOUNDS: {0}--{1}".format(site1, site2))
                    else:
                        logging.info("==> NO COVERAGE: {0}--{1}".format(site1, site2))
                    should_continue = False
        logging.info("==> NUMBER OF COMPLETE REGIONS FOR STARTING SITE {0}: {1}".format(site1, len(complete_regions_for_site1)))
        if len(complete_regions_for_site1) > 0:
            for region in complete_regions_for_site1:
                if(seq_strand == '+'):
                    s = "\t".join([ seq_name, "Cutfiller", "exon",
                                     str(region[ 0 ] + seq_start), str(region[ 1 ] + seq_start), str(region[ 2 ]),
                                     seq_strand, ".",
                                     "exon_id=Ex_Com_" + str(exon_id) + ";"
                                 ])
                else:
                    s = "\t".join([ seq_name, "Cutfiller", "exon",
                                     str(seq_end - region[ 1 ] + 1), str(seq_end - region[ 0 ] + 1), str(region[ 2 ]),
                                     seq_strand, ".",
                                     "exon_id=Ex_Com_" + str(exon_id) + ";"
                                 ])
                print s
                exon_id += 1
    
    logging.info("==> NUMBER OF 3' SPLICE SITES IN WHICH NO COMPLETE REGION STARTS: {0}".format(len(three_p_sites_no_exons)))
    for pos in three_p_sites_no_exons:
        if(cov_reads[ pos - min_aln_pos - args.gapEnds : pos - min_aln_pos + args.gapEnds ].sum() > 0):
            three_p_reg = cov_reads[ pos - min_aln_pos : max_aln_pos - min_aln_pos ]
            three_p_exon = grow_right(three_p_reg, pos, args.gapEnds)
            if(len(three_p_exon) > 0):
                #print "{0}\t{1}\t2".format(three_p_exon[ 0 ], three_p_exon[ 1 ])
                if(seq_strand == '+'):
                    s = "\t".join([ seq_name, "Cutfiller", "exon",
                                     str(pos + seq_start), str(three_p_exon[ 1 ] + seq_start), str(three_p_exon[ 2 ]),
                                     seq_strand, ".",
                                     "exon_id=Ex_3_" + str(exon_id) + ";"
                                 ])
                else:
                    s = "\t".join([ seq_name, "Cutfiller", "exon",
                                     str(seq_end - three_p_exon[ 1 ]), str(seq_end - pos + 1), str(three_p_exon[ 2 ]),
                                     seq_strand, ".",
                                     "exon_id=Ex_3_" + str(exon_id) + ";"
                                 ])
                print s
                exon_id += 1

    logging.info("==> NUMBER OF 5' SPLICE SITES IN WHICH NO COMPLETE REGION ENDS: {0}".format(len(five_p_sites_no_exons)))
    for pos in five_p_sites_no_exons:
        if(cov_reads[ 0 : pos - min_aln_pos + args.gapEnds ].sum() > 0):
            five_p_reg = cov_reads[ 0 : pos - min_aln_pos ]
            five_p_exon = grow_left(five_p_reg, min_aln_pos, args.gapEnds)
            if(len(five_p_exon) > 0):
                #print "{0}\t{1}\t1".format(five_p_exon[ 0 ], five_p_exon[ 1 ])
                if(seq_strand == '+'):
                    s = "\t".join([ seq_name, "Cutfiller", "exon",
                                     str(five_p_exon[ 0 ] + seq_start), str(pos + seq_start), str(five_p_exon[ 2 ]),
                                     seq_strand, ".",
                                     "exon_id=Ex_5_" + str(exon_id) + ";"
                                 ])
                else:
                    s = "\t".join([ seq_name, "Cutfiller", "exon",
                                     str(seq_end - pos + 1), str(seq_end - five_p_exon[ 0 ]), str(five_p_exon[ 2 ]),
                                     seq_strand, ".",
                                     "exon_id=Ex_5_" + str(exon_id) + ";"
                                 ])
                print s
                exon_id += 1

    logging.info("==> PROCESS COMPLETE")
        

if __name__ == "__main__":
    main()

#!/usr/bin/env python2.7

import sys
import argparse
import rbtree
import pysam

def is_unique( interval_set ):
    """
    Analyze interval_set and return true if it represent a single, contiguous
    region.
    """
    uniq = True
    if len( interval_set ) == 0:
        return False
    elif len( interval_set ) == 1:
        print "==> COVERAGE: {0}--{1}".format(interval_set[0].root[0], interval_set[0].root[1])
        return True
    else:
        s = "" + str(interval_set[0].root[0])
        for i in range( len( interval_set ) -1 ):
            if( interval_set[i].root[1] +1 != interval_set[i+1].root[0] ):
                s = s + "--" + str(interval_set[i].root[1]) + "//" + str(interval_set[i+1].root[0])
                uniq = False
        s = s + "--" + str(interval_set[-1].root[1])
        print "==> COVERAGE: {0}".format(s)
        return uniq

def grow_right( interval_set ):
    """
    Compute the biggest contiguous region that starts from the first element of
    interval_set (if any).
    """
    if len( interval_set ) == 0:
        return [ ]
    elif len( interval_set ) == 1:
        return interval_set
    else:
        grown_interval = [ interval_set[ 0 ] ]
        for interval in interval_set[ 1 : ]:
            if interval.root[ 0 ] == grown_interval[ -1 ].root[ 1 ] +1:
                grown_interval.append( interval )
            else:
                return grown_interval
        return grown_interval

def grow_left( interval_set ):
    """
    Compute the biggest contiguous region that ends in the last element of
    interval_set (if any).
    """
    if len( interval_set ) == 0:
        return [ ]
    elif len( interval_set ) == 1:
        return interval_set
    else:
        grown_interval = [ interval_set[ -1 ] ]
        for interval in interval_set[ 0:-1 ][::-1]:
            if grown_interval[ -1 ].root[ 0 ] == interval.root[ 1 ] +1 :
                grown_interval.append( interval )
            else:
                return grown_interval[::-1]
        return grown_interval[::-1]

def main( ):
    parser = argparse.ArgumentParser( prog = "cutfiller",
                                      description = "Rebuild exons from PIntron output.",
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument( '-p', '--pintron-output', help = "PIntron output file",
                         required = True, dest = 'pfile' )
    parser.add_argument( '-a', '--alignment-file', help = 'Alignment file (PIntron)',
                         required = False, dest = 'alfile')
    parser.add_argument( '-s', '--samfile', help = 'Alignment file (sam format)',
                         required = False, dest = 'samfile')
    parser.add_argument( '-l', '--max-intron-length',
                         help = 'Max intron length accepted. Introns that exceed this value will be discarded',
                         required = False, dest = 'maxIntron', type = int, default = 15000 )
    parser.add_argument( '-g', '--gap-ends', help = 'Gap threshold at exon ends for their computation.',
                         required = False, dest = 'gapEnds', type = int, default = 3)
    args = parser.parse_args(  )

    if not args.alfile and not args.samfile:
        print >> sys.stderr, "ERROR: No samfile given and no PIntron alignment file given.\nAborting..."
        sys.exit( 1 )
    if args.alfile and args.samfile:
        print >> sys.stderr, "ERROR: Both samfile and PIntron alignment files given.\nAborting..."
        sys.exit( 1 )

    cutst = rbtree.RBIntervalTree( )

    print "==> OPENING ALIGNMENT FILE"

    # Open PIntron alignment file and build the tree
    align_index = 0
    min_aln_pos = 0
    max_aln_pos = 0
    if args.alfile:
        with open( args.alfile, 'r' ) as inInts:
            for line in inInts.readlines( ):
                if not line.startswith( ">" ) and not line.startswith( "#" ):
                    print >> sys.stderr, "{0:<50}\r".format( "==> PROCESSING ALIGNMENT NUMBER {0:<10}".format( align_index ) ),
                    elements = line.split( ' ' )
                    begin = int( elements[ 2 ] )
                    end = int( elements[ 3 ] )
                    if(align_index == 0):
                        min_aln_pos = begin
                        max_aln_pos = end
                    if begin < min_aln_pos:
                        min_aln_pos = begin
                    if end > max_aln_pos:
                        max_aln_pos = end
                    align_index += 1
                    cutst.rbinsert( [ begin, end ] )
    elif args.samfile:
        with open( args.samfile, 'r') as inSam:
            samfile = pysam.Samfile( args.samfile )
            for align in samfile:
                print >> sys.stderr, "{0:<50}\r".format( "==> PROCESSING ALIGNMENT NUMBER {0:<10}".format( align_index ) ),
                align_index += 1
                begin = int( align.pos )
                end = begin + int( align.qlen )
                cutst.rbinsert( [begin, end] )

    print >> sys.stderr # newline

    five_p_sites = []
    five_p_sites_no_exons = []
    three_p_sites = []
    three_p_sites_no_exons = []
    
    print "==> OPENING SPLICE SITE FILE: {0}".format( args.pfile)
    with open( args.pfile, 'r' ) as inIntrons:
        for line in inIntrons.readlines( ):
            elements = line.split( "\t" )
            region_end = int( elements[ 0 ] ) -1
            region_begin = int( elements[ 1 ] ) +1
            if abs( region_begin - region_end ) <  args.maxIntron:
                if region_end not in five_p_sites:
                    five_p_sites.append( region_end )
                    five_p_sites_no_exons.append ( region_end )
                if region_begin not in three_p_sites:
                    three_p_sites.append( region_begin )
                    three_p_sites_no_exons.append( region_begin )
            else:
                print "==> INTRON TOO LARGE: [{0}, {1}]".format( region_begin, region_end )

    three_p_sites = sorted( list( set( three_p_sites ) ) )
    five_p_sites = sorted( list( set( five_p_sites ) ) )

    five_p_sites_no_exons = sorted( list( set( five_p_sites_no_exons ) ) )
    three_p_sites_no_exons = sorted( list( set( three_p_sites_no_exons ) ) )

    print "==> PROCESSING 3' SPLICE SITES"
    
    # Analyze complete regions for every splice site
    for site1 in three_p_sites:
        print "==> COMPUTING NUMBER OF COMPLETE REGIONS FOR STARTING SITE {0}".format( site1 )
        complete_regions_for_site1 = []
        should_continue = True
        for site2 in five_p_sites:
            if should_continue and site1 < site2:
                nodes = cutst.search( [ site1, site2 ] )
                if ( len( nodes ) > 0 and
                     abs( nodes[ 0 ].root[ 0 ] - site1 ) < args.gapEnds and
                     abs( nodes[ -1 ].root[ 1 ] - site2 ) < args.gapEnds
                ):
                    print "==> INTERVAL: {0}--{1}".format(site1, site2)
                    if ( is_unique( nodes ) ):
                        if site1 in three_p_sites_no_exons:
                            three_p_sites_no_exons = filter( lambda x: x != site1, three_p_sites_no_exons )
                        if site2 in five_p_sites_no_exons:
                            five_p_sites_no_exons = filter( lambda x: x != site2, five_p_sites_no_exons )
                        complete_regions_for_site1.append( (site1, site2) )
                    else:
                        should_continue = False
                else:
                    print "==> INTERVAL: {0}--{1}".format(site1, site2)
                    if ( len( nodes ) > 0):
                        print "==> EXCEEDING BOUNDS"
                    else:
                        print "==> NO COVERAGE"
                    should_continue = False
        print "==> NUMBER OF COMPLETE REGIONS FOR STARTING SITE {0}: {1}".format( site1, len(complete_regions_for_site1)  )
        if len( complete_regions_for_site1 ) > 0:
            for region in complete_regions_for_site1:
                print "{0}\t{1}\t0".format(*region)

    print "==> NUMBER OF 3' SPLICE SITES IN WHICH NO COMPLETE REGION STARTS: {0}".format( len( three_p_sites_no_exons ) )

    for position in three_p_sites_no_exons:
        nodes = cutst.search( [ position, cutst.get_max( ) ] )
        nodes = grow_right( nodes )
        if len( nodes ) >= 1:
            print "{0}\t{1}\t2".format( position, nodes[ -1 ].root[ 1 ] )

    print "==> NUMBER OF 5' SPLICE SITES IN WHICH NO COMPLETE REGION ENDS: {0}".format( len( five_p_sites_no_exons ) )

    for position in five_p_sites_no_exons:
        nodes = cutst.search( [ 0, position ] )
        nodes = grow_left( nodes )
        if len( nodes ) >= 1:
            print "{0}\t{1}\t1".format( nodes[ 0 ].root[ 0 ], position )
        nodes = [ ]

    print "==> PROCESS COMPLETE"
    
if __name__ == "__main__":
    main( )

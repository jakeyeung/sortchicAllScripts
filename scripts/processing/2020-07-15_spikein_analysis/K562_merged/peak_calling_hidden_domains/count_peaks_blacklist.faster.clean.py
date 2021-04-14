#!/usr/bin/env python
'''
DESCRIPTION

    Counts peaks nonpeaks blacklist by chromosome

FOR HELP

    python count_peaks_nonpeaks_blacklist.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-11-07
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import csv
import sys, argparse, datetime
import singlecellmultiomics.modularDemultiplexer
import os
import sys
import pysam
import collections
import argparse
import pandas as pd
import numpy as np
import itertools
import singlecellmultiomics.modularDemultiplexer.baseDemultiplexMethods
import gzip  # for loading blacklist bedfiles

def read_has_alternative_hits_to_non_alts(read):
    if read.has_tag('XA'):
        for alt_align in read.get_tag('XA').split(';'):
            if len(alt_align) == 0:  # Sometimes this tag is empty for some reason
                continue

            hchrom, hpos, hcigar, hflag = alt_align.split(',')
            if not hchrom.endswith('_alt'):
                return True
    return False

def read_should_be_counted(read, args, verbose=False):
    """
    Check if a read should be counted given the filter arguments

    Parameters
    ----------
    read : pysam.AlignedSegment or None
        read to check if it should be counted

    Returns
    -------
    bool
    """

    if args.r1only and read.is_read2:
        if verbose:
            print("read2 but r1only")
        return False
    if args.r2only and read.is_read1:
        if verbose:
            print("read1 but r2only")
        return False

    if args.filterMP:
        if not read.has_tag('mp'):
            if verbose:
                print("Not mp")
            return False
        if read.get_tag('mp')=='unique':
            return True
        return False

    if read is None or read.is_qcfail:
        if verbose:
            print("QC fail")
        return False

    # Mapping quality below threshold
    if read.mapping_quality < args.minMQ:
        if verbose:
            print("Less than minMQ", args.minMQ)
        return False


    if args.proper_pairs_only and not read.is_proper_pair:
        if verbose:
            print("not proper pair")
        return False

    if args.no_indels and ('I' in read.cigarstring or 'D' in read.cigarstring):
        if verbose:
            print("Has indels")
        return False

    if args.max_base_edits is not None and read.has_tag('NM') and int(read.get_tag('NM'))>args.max_base_edits:
        if verbose:
            print("Too many base edit")
        return False

    if args.no_softclips and 'S' in read.cigarstring:
        if verbose:
            print("Has softclips")
        return False

    # Read has alternative hits
    if args.filterXA:
        if read_has_alternative_hits_to_non_alts(read):
            if verbose:
                print("Has alt hits")
            return False

    # Read is a duplicate
    # (args.dedup and ( not read.has_tag('RC') or (read.has_tag('RC') and read.get_tag('RC')!=1))):
    if read.is_unmapped or \
        (args.dedup and (read.has_tag("RR") or read.is_duplicate)):
        if verbose:
            return False
    return True


def read_in_dic(read, dic, debug=False):
    if read.reference_name in dic:
        if debug:
            print("Read ref name in dic:", read.reference_name)
        # iterate through tuples of startend to check
        for startend in dic[read.reference_name]:
            if debug:
                print("Startend:", startend)
                print("Read ref start:", read.reference_start)
                print("Read ref end", read.reference_end)
            # start is 0-based inclusive, end is 0-based exclusive
            start_inside = read.reference_start >= startend[0] and read.reference_start < startend[1]
            end_inside = read.reference_end >= startend[0] and read.reference_end < startend[1]
            # start_bad = read.reference_start >= startend[0] and read.reference_start < startend[1]
            # end_bad = read.reference_end >= startend[0] and read.reference_end < startend[1]
            if debug:
                print("Start inside:", start_inside)
                print("End inside:", end_inside)
                print("Start or End inside:", start_inside or end_inside)
                # input("Waiting...")
            if start_inside or end_inside:
                return True
            # else:
            #     return False
    else:
        # read was not in peak
        return False


def index_bedfile(infbed):
    # create blacklist dictionary {chromosome : [ (start1, end1), ..., (startN, endN) ]}
    # used to check each read and exclude if it is within any of these start end sites
    #
    bedlist_dic = {}
    print("Creating list dictionary:")
    with open(infbed, mode='r') as blfile:
        for row in blfile:
            parts = row.strip().split()
            chromo, start, end = parts[0], int(parts[1]), int(parts[2])
            if chromo not in bedlist_dic:
                # init chromo
                bedlist_dic[chromo] = []  # list of tuples
            bedlist_dic[chromo].append( (start, end) )
    return(bedlist_dic)
    # print(bedlist_dic)

def main():
    parser = argparse.ArgumentParser(description='Counts peaks nonpeaks blacklist by chromosome')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam file')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output samples, counts, chromosome')
    parser.add_argument('-bedfile', metavar='bedfile of peaks',
                        help='Bed of peaks')
    parser.add_argument('-mode', metavar='cuts_in_peak or cuts_in_chromo',
                        help='cuts_in_peak or cuts_in_chromo')
    parser.add_argument('-blacklist', metavar='bedfile of blacklist',
                        help='Blacklist')
    parser.add_argument('-minMQ', metavar='INT', type=int, default=0,
                        help='Minimum mapping quality')
    parser.add_argument('-max_base_edits', metavar='INT', type=int,
                        help='Count reads with at most this value of bases being different than the reference')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--dedup', action='store_true',
                        help='Count only the first occurence of a molecule. Requires RC tag to be set. Reads without RC tag will be ignored!')
    parser.add_argument('--proper_pairs_only', action='store_true',
                        help='Only count reads mapped in a proper pair (within the expected insert size range)')
    parser.add_argument('--no_softclips', action='store_true',
                        help='Only count reads without softclips')
    parser.add_argument('--no_indels', action='store_true',
                        help='Only count reads without indels')
    parser.add_argument('--r1only', action='store_true',
                        help='Count r1only')
    parser.add_argument('--r2only', action='store_true',
                        help='Count r2only')
    parser.add_argument('--filterXA', action='store_true',
                        help='Filter XA (alternative hits) tag')
    parser.add_argument('--filterMP', action='store_true',
                        help='Filter reads which are not uniquely mappable, this is based on presence on the `mp` tag, which can be generated by bamtagmultiome.py')
    parser.add_argument('--debug', action='store_true',
                        help='Debug mode, default False')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    debug = args.debug

    # set up count dic
    count_dic = {}

    # set up peak dic
    if args.mode == "cuts_in_peak":
        print("Indexing peak list")
        peaks_dic = index_bedfile(args.bedfile)

    # create blacklist dictionary {chromosome : [ (start1, end1), ..., (startN, endN) ]}
    # used to check each read and exclude if it is within any of these start end sites
    #
    blacklist_dic = index_bedfile(args.blacklist)

    assert args.mode == "cuts_in_peak" or args.mode == "cuts_in_chromo"
        
    with pysam.AlignmentFile(args.infile, "rb") as samfile:
        if args.mode == "cuts_in_peak":
            for chromo, startend_tuplst in peaks_dic.items():
                for startend in startend_tuplst:
                    if debug:
                        print('startend:', startend)
                        print('startend0:', startend[0])
                        print('startend1:', startend[1])
                        input("Waiting....")
                    for read in samfile.fetch(contig = chromo, start = startend[0], stop = startend[1]):
                        samp = read.get_tag('SM')
                        readchromo = read.reference_name
                        if samp not in count_dic:
                            # count_dic[samp] = {'cuts_in_peak' : 0, 'cuts_in_blacklist' : 0, 'skipped' : 0}
                            count_dic[samp] = {}
                        if readchromo not in count_dic[samp]:
                            count_dic[samp][readchromo] = {'cuts_in_peak' : 0, 'cuts_in_blacklist' : 0, 'skipped' : 0}

                        if not read_should_be_counted(read, args):
                            count_dic[samp][readchromo]['skipped'] += 1 
                            continue
                        # Read is in blacklist
                        if debug:
                            print("Checking blacklist")
                        if read_in_dic(read, blacklist_dic, debug = debug):
                            count_dic[samp][readchromo]['cuts_in_blacklist'] += 1
                            continue
                        count_dic[samp][readchromo]['cuts_in_peak'] += 1
                        if debug:
                            print(count_dic[samp][readchromo])
                            input("Waiting2...")
                            # skip and write to output

        elif args.mode == "cuts_in_chromo":
            for read in samfile.fetch():
                samp = read.get_tag('SM')
                readchromo = read.reference_name
                if samp not in count_dic:
                    # count_dic[samp] = {'cuts_in_chromo' : 0, 'cuts_in_blacklist' : 0, 'skipped' : 0}
                    count_dic[samp] = {}
                if readchromo not in count_dic[samp]:
                    count_dic[samp][readchromo] = {'cuts_in_chromo' : 0, 'cuts_in_blacklist' : 0, 'skipped' : 0}

                if not read_should_be_counted(read, args):
                    count_dic[samp][readchromo]['skipped'] += 1 
                    continue
                # Read is in blacklist
                if debug:
                    print("Checking blacklist")
                if read_in_dic(read, blacklist_dic, debug = debug):
                    count_dic[samp][readchromo]['cuts_in_blacklist'] += 1
                    continue
                count_dic[samp][readchromo]['cuts_in_chromo'] += 1
        else:
            print("Warning: args.mode should be cuts_in_peak or cuts_in_chromo", args.mode)
    # write count_dic as table
    with open(args.outfile, "w") as outf:
        jwriter = csv.writer(outf, delimiter = "\t")
        # write header
        if args.mode == "cuts_in_peak":
            jheader = ['samp', 'chromosome', 'cuts_in_peak', 'cuts_in_blacklist', 'bad_reads']
        elif args.mode == "cuts_in_chromo":
            jheader = ['samp', 'chromosome', 'cuts_in_chromo', 'cuts_in_blacklist', 'bad_reads']
        else:
            print("Warning: args.mode should be cuts_in_peak or cuts_in_chromo")
        jwriter.writerow(jheader)
        for samp, chromodic in count_dic.items():
            for chromo, count_subdic in chromodic.items():
                if args.mode == "cuts_in_peak":
                    outrow = [samp, chromo, count_subdic['cuts_in_peak'], count_subdic['cuts_in_blacklist'], count_subdic['skipped']]
                elif args.mode == "cuts_in_chromo":
                    outrow = [samp, chromo, count_subdic['cuts_in_chromo'], count_subdic['cuts_in_blacklist'], count_subdic['skipped']]
                else:
                    print(args.mode)
                    print("Warning: argsmode should be cuts_in_peak or cuts_in_chromo", args.mode)
                jwriter.writerow(outrow)

if __name__ == '__main__':
    main()

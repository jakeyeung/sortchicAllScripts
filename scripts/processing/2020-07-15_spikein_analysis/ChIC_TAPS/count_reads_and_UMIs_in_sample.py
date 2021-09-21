#!/usr/bin/env python
'''
DESCRIPTION

    Get UMI and reads by sample

FOR HELP

    python count_reads_and_UMIs_in_sample.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-09-25
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import csv

def main():
    parser = argparse.ArgumentParser(description='Get UMI and reads by sample')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output text file')
    parser.add_argument('-mapq', metavar='MAPQ', type=int, default=40, 
                        help='Minimum MAPQ ')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
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

    samp_dic = {}
    with pysam.AlignmentFile(args.infile, "rb") as inbam:
        for read in inbam.fetch():
            if read.mapping_quality < args.mapq:
                continue
            samp = read.get_tag("SM")
            # init if needed
            if samp not in samp_dic:
                samp_dic[samp] = {"read" : 0, "umi" : 0}  # reads and umi tuple
            samp_dic[samp]["read"] += 1
            if not read.is_duplicate:
                samp_dic[samp]["umi"] += 1
    # write to output
    with open(args.outfile, "w") as outf:
        outwriter = csv.writer(outf, delimiter = "\t")
        for samp, rdic in samp_dic.items():
            outrow = [samp, rdic["read"], rdic["umi"]]
            outwriter.writerow(outrow)

if __name__ == '__main__':
    main()

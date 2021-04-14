#!/usr/bin/env python
'''
DESCRIPTION

    Filter for only one chromosome

FOR HELP

    python split_bam_by_chromosome.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-09-20
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import os
import pysam

def main():
    parser = argparse.ArgumentParser(description='Filter for only one chromosome')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output bam')
    parser.add_argument('-chromo', metavar='CHROMOSOME',
                        help='Chromosome to filter out')
    parser.add_argument('--dedup', action='store_true', 
                        help='Deduplicate reads')
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

    assert not os.path.exists(args.outfile)
    assert os.path.exists(args.infile)


    with pysam.AlignmentFile(args.infile, "rb") as inbam:
        outbam = pysam.AlignmentFile(args.outfile, "wb", template = inbam)
        for read in inbam.fetch(args.chromo):
            if read.is_duplicate and args.dedup:
                continue
            outbam.write(read)
        outbam.close()    
    pysam.index(args.outfile)


if __name__ == '__main__':
    main()

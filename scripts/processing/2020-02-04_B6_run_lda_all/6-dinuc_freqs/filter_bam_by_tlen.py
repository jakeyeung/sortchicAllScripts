#!/usr/bin/env python
'''
DESCRIPTION

    Filter bam by tlen to get precise dyads

FOR HELP

    python filter_bam_by_tlen.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-08-05
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import os
import sys, argparse, datetime
import pysam

def main():
    parser = argparse.ArgumentParser(description='Filter bam by tlen to get precise dyads')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output bam')
    parser.add_argument('-length_min', metavar='template length', type=int,
                        help='minimum length')
    parser.add_argument('-length_max', metavar='template length', type=int,
                        help='maximum length')
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

    inbam = pysam.AlignmentFile(args.infile, "rb")
    header = inbam.header.copy()
    outbam = pysam.AlignmentFile(args.outfile, "wb", header=header)

    nkept = 0
    ntotal = 0

    for read in inbam.fetch():
        ntotal += 1
        tlen = read.template_length
        if tlen >= args.length_min and tlen <= args.length_max:
            nkept += 1
            outbam.write(read)


    outbam.close()
    inbam.close()

    cmd="samtools index %s" % args.outfile
    os.system(cmd)

    print("Readskept / Readstotal")
    print(nkept, "/", ntotal)


if __name__ == '__main__':
    main()

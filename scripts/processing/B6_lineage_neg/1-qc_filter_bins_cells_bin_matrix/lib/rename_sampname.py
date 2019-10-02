#!/usr/bin/env python
'''
DESCRIPTION

    Description

FOR HELP

    python rename_sampname.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2019-05-10
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam


def main():
    parser = argparse.ArgumentParser(description='Description')
    parser.add_argument('infile', metavar='INFILE',
                        help='bam input')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='bam output')
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

    with pysam.AlignmentFile(args.infile, "rb") as bamf:
        with pysam.AlignmentFile(args.outfile, "wb", template=bamf) as outbam:
            for read in bamf:
                # Is:NS500813;RN:490;Fc:H5W5KBGXB;La:1;Ti:12209;CX:23210;CY:11779;Fi:N;CN:0;aa:GTGAAA;aA:GTGAAA;aI:19;LY:B6-13W1-BM-H3K4me1-1_AH5W5KBGXB_S1;RX:CAG;RQ:GGG;BI:160;bc:TCTATCTC;BC:TCTATCTC;QT:GGKKKKKK;MX:NLAIII384C8U3
                headerlst = read.query_name.split(";")
                sampname = headerlst[12]
                # clip _AH5W5KBGXB_S1
                sampnamenew = sampname.split("_")[0]
                headerlst[12] = sampnamenew
                querynamenew = ';'.join(headerlst)
                read.query_name = querynamenew
                outbam.write(read)
    # sort bam file
    pysam.index(args.outfile)


if __name__ == '__main__':
    main()

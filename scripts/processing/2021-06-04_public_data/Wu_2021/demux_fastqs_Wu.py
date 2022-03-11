#!/usr/bin/env python
'''
DESCRIPTION

    Read r1 and r2, move barcode into r2

FOR HELP

    python demux_fastqs_Ku.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2021-06-14
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import csv
import gzip


def main():
    parser = argparse.ArgumentParser(description='Read r1 and r2, move barcode into r1')
    parser.add_argument('-infile1', metavar='FASTQ FILE 1',
                        help='fastq1 r1')
    parser.add_argument('-infile2', metavar='FASTQ FILE 2',
                        help='fastq2 r2')
    parser.add_argument('-outfile', metavar='OUTPUT FASTQ FILE 2',
                        help='Output fastq2 r1 info appended')
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

    # read r1
    fq1 = pysam.FastxFile(args.infile1, "r")
    fq2 = pysam.FastxFile(args.infile2, "r")
    fq2out = gzip.open(args.outfile, mode="wt")

    ntotal = 0
    r1length = 25
    splitindx = 11  # first 11 bp are inline, next 14 bp are UMIs

    with gzip.open(args.outfile, mode = "wt") as fq2out:
        for f1, f2 in zip(fq1, fq2):
            ntotal += 1
            # https://auto-process-ngs.readthedocs.io/en/latest/single_cell/icell8.html#icell8-single-cell-atac-seq
            # first 11 bp are inline barcode, the next 14 bp are UMIs
            bc_umi = f1.sequence
            assert len(bc_umi) == r1length

            bc = f1.sequence[0:splitindx]
            umi = f1.sequence[splitindx:r1length]

            bc_umi_toappend = '_'.join([bc, umi])

            f2new = f2
            f2new.name = '.'.join([f2new.name, bc_umi_toappend])

            # write new f2 to file 
            fq2out.write(str(f2new))
            fq2out.write("\n")
    fq1.close()
    fq2.close()

    print("ntotal=", ntotal)

if __name__ == '__main__':
    main()

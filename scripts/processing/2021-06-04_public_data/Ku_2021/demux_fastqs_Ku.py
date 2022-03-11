#!/usr/bin/env python
'''
DESCRIPTION

    Read r1 and r2, move barcode into r1

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
    parser.add_argument('-bcfile', metavar='BARCODES TXT',
                        help='Each row is cell barcode')
    parser.add_argument('-outfile', metavar='OUTPUT FASTQ FILE 1',
                        help='Output fastq1 r1 appended')
    parser.add_argument('-bclength', type = int, metavar='INT', default=8,
                        help='integer')
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

    # read barcodes
    bclength = 8
    # bcs = []
    bcs = {}
    with open(args.bcfile, "r") as bcf:
        bcfreader = csv.reader(bcf, delimiter = "\t")
        for row in bcfreader:
            bc = row[0]
            assert len(bc) == bclength
            bcs[bc] = True
            # bcs.append(bc)

    print(len(bcs), " barcodes read into memory")
    print(bcs)
    bcs_tup = tuple(bcs)

    # read r1
    fq1 = pysam.FastxFile(args.infile1, "r")
    fq2 = pysam.FastxFile(args.infile2, "r")
    fq1out = gzip.open(args.outfile, mode="wt")

    nbad = 0
    ntotal = 0

    with gzip.open(args.outfile, mode = "wt") as fq1out:
        for f1, f2 in zip(fq1, fq2):
            ntotal += 1
            f2bc = f2.sequence[:bclength]
            if f2bc in bcs:
                f1new = f1
                f1new.name = '.'.join([f1.name, f2bc])
                # write new f1 to file 
                fq1out.write(str(f1new))
                fq1out.write("\n")
                # print(f1.name)
                # print(f1new.name)
                # print(f2.sequence)
                # print(f2bc)
                # print("Written")
            else:
                nbad += 1
                # print(f2.sequence)
                # print("Skipped")
            # input()

    fq1.close()
    fq2.close()

    print("ntotal=", ntotal)
    print("nbad=", nbad)

    print("Frac of reads thrown:", float(nbad) / float(ntotal))

    # bashCommand="".join(["gzip ", args.outfile])
    # os.system(bashCommand)


if __name__ == '__main__':
    main()

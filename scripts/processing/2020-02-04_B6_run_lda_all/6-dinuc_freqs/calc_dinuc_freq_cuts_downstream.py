#!/usr/bin/env python
'''
DESCRIPTION

    Given a list of sequences, calculate dinucleotide frequency

FOR HELP

    python calc_dinuc_freq_cuts_downstream.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-08-05
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import os
import sys, argparse, datetime
import gzip
import csv

def main():
    parser = argparse.ArgumentParser(description='Given a list of sequences, calculate dinucleotide frequency')
    parser.add_argument('infile', metavar='INFILE',
                        help='Tab csv (zipped), second column is sequence')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output dinuc frequency')
    parser.add_argument('-baseseq', metavar='sequence', default = 'AT',
                        help='Nucloetide frequency to check. Default AT')
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

    with open(args.outfile, mode="w") as outf:
        outwriter = csv.writer(outf, delimiter = "\t")
        if args.infile.endswith(".gz"):
            f = gzip.open(args.infile, mode="rt")
        else:
            f = open(args.infile, mode="r")
        
        freader = csv.reader(f, delimiter = "\t")
        for row in freader:
            seq = row[1]
            # iterate sequence to calculate dinuc freq
            dinuc = []
            for i in range(len(seq) - 1):
                subseq = seq[i : (i + 2)]
                if subseq[0] in args.baseseq and subseq[1] in args.baseseq:
                    dinuc.append(1)
                else:
                    dinuc.append(0)
            outwriter.writerow(dinuc)
        f.close()
    cmd=" ".join(["gzip", args.outfile])
    os.system(cmd)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
DESCRIPTION

    Input Tss and Tes, output bed file

FOR HELP

    python create_tss_tes_bedfile.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-08-10
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import gzip
import csv

def main():
    parser = argparse.ArgumentParser(description='Input Tss and Tes, output bed file')
    parser.add_argument('-TssBed', metavar='BEDFILE',
                        help='TSS file')
    parser.add_argument('-TesBed', metavar='BEDFILE',
                        help='TES file')
    parser.add_argument('-OutBed', metavar='BEDFILE',
                        help='Output bed file')
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

    # load annots
    annots = {}
    with gzip.open(args.TssBed, "rt") as tssf:
        jreader = csv.reader(tssf, delimiter = "\t")
        for row in jreader:
            jchromo, jstart, jend, strand, jname = row[0], row[1], row[2], row[3], row[4]
            if jname not in annots:
                annots[jname] = [jchromo, strand, jstart]  # add TES later
            else:
                print("TSS Skipping duplicated:", jname)
                continue

    with gzip.open(args.TesBed, "rt") as tesf:
        jreader = csv.reader(tesf, delimiter = "\t")
        for row in jreader:
            jchromo, jstart, jend, strand, jname = row[0], row[1], row[2], row[3], row[4]
            if jname not in annots:
                print("TES name not found, skipping:", jname)
                continue
            else:
                annots[jname].append(jstart)

    # write output 
    with gzip.open(args.OutBed, "wt") as outf:
        jwriter = csv.writer(outf, delimiter = "\t")
        for jname, jlist in annots.items():
            if len(jlist) != 4:
                print("TSS and TES incomplete, skipping:", jname)
                continue
            jchromo, strand, tss, tes = jlist[0], jlist[1], jlist[2], jlist[3]
            if strand == "+":
                jstart, jend = tss, tes
            else:
                jstart, jend = tes, tss
            assert int(jend) - int(jstart) > 0
            outrow = [jchromo, jstart, jend, strand, jname]
            jwriter.writerow(outrow)
        

if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
DESCRIPTION

    Remove invalid records after HD

FOR HELP

    python remove_invalid_records.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2021-06-06
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import csv
import gzip

def main():
    parser = argparse.ArgumentParser(description='Remove invalid records after HD')
    parser.add_argument('infile', metavar='INFILE',
                        help='hiddendomains analysis bed')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='hiddendomains filtered analysis bed')
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

    badcount = 1

    with open(args.infile, "r") as infile:
        jreader = csv.reader(infile, delimiter = "\t")
        with open(args.outfile, "w") as outfile:
            jwriter = csv.writer(outfile, delimiter = "\t")
            for row in jreader:
                jstart = int(row[1])
                jend = int(row[2])
                if jend < jstart: 
                    badcount += 1
                else:
                    outrow = row
                    jwriter.writerow(outrow)
    print("Removed", badcount, "rows")

if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
DESCRIPTION

    Count reads from a bedfile before and after filtering for ones that fall in peaks

FOR HELP

    python count_reads_from_bed_before_and_after_filtering.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2021-06-06
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import csv
import gzip

def openfile(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)



def main():
    parser = argparse.ArgumentParser(description='Count reads from a bed fie ')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bed before filtering. Can be gz')
    parser.add_argument('infilefilt', metavar='INFILE_FILTERED',
                        help='Input bed after filtering')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output cell name, counts in file 1, counts in file 2')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--writeheader', action='store_true',
                        help='Write header otherwise no header')
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

    with openfile(args.infile, 'r') as infile:
        for i, line in enumerate(infile):
            # i starts at 0
            pass
        nlines1 = i + 1

    with openfile(args.infilefilt, 'r') as infilefilt:
        for i2, line2 in enumerate(infilefilt):
            # i starts at 0
            pass
        nlines2 = i2 + 1

    # save output
    jheader = ['infile1', 'infile2', 'nlines1', 'nlines2']
    outrow = [args.infile, args.infilefilt, nlines1, nlines2]

    with open(args.outfile, 'w') as outfile:
        outobj = csv.writer(outfile, delimiter = "\t")
        if args.writeheader:
            outobj.writerow(jheader)
        outobj.writerow(outrow)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
DESCRIPTION

    Run bamCoverage with custom scale factor: multiplies scale factor by TotalCounts so it cancels out the default normalization factor

FOR HELP

    python calculate_scale_factors_for_deeptools.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-06-19
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
from deeptools import utilities, bamHandler

def main():
    parser = argparse.ArgumentParser(description='Run bamCoverage with custom scale factor')
    parser.add_argument('-infile', metavar='INFILE',
                        help='Input bam file for bamCoverage')
    parser.add_argument('-outfile', metavar='OUTFILE',
                        help='Output bigwig file')
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

    # calculate the default normalization factor
    bam_handle, mapped, unmapped, stats = bamHandler.openBam(args.infile, returnStats=True)

    print(bam_handle)

    bam_mapped_total = utilities.bam_total_reads(bam_handle, False, stats)

    print(bam_mapped_total)

    with open(args.outfile, "w") as outf:
        outf.write(str(bam_mapped_total))

if __name__ == '__main__':
    main()

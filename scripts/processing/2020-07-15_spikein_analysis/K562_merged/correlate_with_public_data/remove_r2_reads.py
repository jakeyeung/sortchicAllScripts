#!/usr/bin/env python
'''
DESCRIPTION

    Remove r2 reads

FOR HELP

    python remove_r2_reads.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-11-27
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''
import os 
import sys, argparse, datetime
import pysam

def main():
    parser = argparse.ArgumentParser(description='Remove r2 reaads')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output bam, sorted and indexed')
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

    outfiletmp=''.join([args.outfile, ".tmpbam"])
    r2count = 0 

    with pysam.AlignmentFile(args.infile, "rb") as infobj:
        with pysam.AlignmentFile(outfiletmp, "wb", template = infobj) as outfobj:
            for read in infobj:
                if read.is_read1:
                    outfobj.write(read)
                else:
                    r2count += 1

    print("Number of r2 reads: %s" % r2count)
    # sort and index
    pysam.sort('-o', args.outfile, outfiletmp)
    pysam.index(args.outfile)

    os.remove(outfiletmp)

if __name__ == '__main__':
    main()

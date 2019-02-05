#!/usr/bin/env python
'''
DESCRIPTION

    Given a list of bams and associated clusters, merge bams by cluster.

FOR HELP

    python merge_bams_by_cluster.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2018-12-23
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import csv, pysam


def main():
    parser = argparse.ArgumentParser(description='Given a list of bams and associated clusters, merge bams by cluster.')
    parser.add_argument('indir', metavar='INDIR',
                        help='Dir of bams')
    parser.add_argument('bam_annot', metavar='CLSTRFILE',
                        help='Path to textfile linking bams to a cluster. Expect full paths in this file in first column, cluster ID in second column.')
    parser.add_argument('outdir', metavar='OUTDIR',
                        help='Dir of bams output')
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

    clstbam = {}
    clstrs = set()
    inbam = {}

    # read bam file
    with open(args.bam_annot, "rb") as rf:
        reader = csv.reader(rf, delimiter = "\t")
        # expects column name
        header = next(reader)
        for row in reader:
            bampath, cluster = row[0], row[1]
            if cluster not in clstbam:
                clstbam[cluster] = []
            clstbam[cluster].append(bampath)
            clstrs.add(cluster)
    clstrs = sorted(list(clstrs))

    # get output files

    # merge bam file
    for c in clstrs:
        outbamfile = pysam.AlignmentFile(outbam[c], "wb")
        with pysam.AlignmentFile(clstbam[c], "rb") as inbamfile:



if __name__ == '__main__':
    main()

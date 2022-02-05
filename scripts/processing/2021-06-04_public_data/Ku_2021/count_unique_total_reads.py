#!/usr/bin/env python
'''
DESCRIPTION

    Count UMIs and total reads for each cell. Output for each cell barcode

FOR HELP

    python count_unique_total_reads.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2021-06-15
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import csv
import gzip


def main():
    parser = argparse.ArgumentParser(description='Count UMIs and total reads for each cell. Output for each cell barcode')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam')
    parser.add_argument('outfile', metavar='OUTFILE GZIP',
                        help='Output unique reads and total reads for each cellbarcode')
    parser.add_argument('-species', metavar='human', default = "human",
                        help='Which chromosomes to keep. Default human')
    parser.add_argument('-mapq', metavar='integer', type=int, default = 40,
                        help='Min mapping quality. Default 40')
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

    assert args.outfile.endswith(".gz")

    # chromos keep 
    chromos = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chromos = set(chromos)


    # track uniq reads and position for every read 
    print("Collecting dic")
    outdic = {}
    skipped = 0
    total = 0
    with pysam.AlignmentFile(args.infile, "rb") as bamfile:
        for entry in bamfile:
            # print(entry.query_name)
            total += 1
            umi = entry.query_name.split(".")[2]
            start = entry.reference_start
            end = entry.reference_end
            chromo = entry.reference_name
            if chromo not in chromos or entry.mapping_quality <= args.mapq:
                skipped += 1
                continue
            # print(chromo, start, end)
            coord = '.'.join([chromo, str(start), str(end)])
            if umi not in outdic:
                outdic[umi] = {}
            if coord not in outdic[umi]:
                outdic[umi][coord] = 1
            else:
                outdic[umi][coord] += 1
            # print(umi, chromo, start, end)
            # input("waiting...")

    print("Total:", total)
    print("Skipped:", skipped)
    print("Frac skipped", float(skipped) / float(total))

    print("Writing dic")
    # write outdic as bedfile 
    # coord, start, end, dupCount, cellname
    with gzip.open(args.outfile, "wt") as outf:
        outwriter = csv.writer(outf, delimiter = "\t")
        for umi, subdic in outdic.items():
            for coord, dupcount in subdic.items():
                chromo = coord.split(".")[0]
                start = coord.split(".")[1]
                end = coord.split(".")[2]
                outrow = [chromo, start, end, dupcount, umi]
                outwriter.writerow(outrow)


if __name__ == '__main__':
    main()

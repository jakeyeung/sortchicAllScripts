#!/usr/bin/env python
'''
DESCRIPTION

    Get location of counts within window

FOR HELP

    python calculate_cut_positions.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-08-06
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import gzip 
import csv

def main():
    parser = argparse.ArgumentParser(description='Get location of counts within window')
    parser.add_argument('-infile', metavar='INFILE',
                        help='Input bam file')
    parser.add_argument('-outprefix', metavar='OUTPREFIX',
                        help='Output twofiles: (1)  a matrix one row one region, (2) number of columns is radius * 2 + 1')
    parser.add_argument('-bedfile', metavar='BEDFILE',
                        help='Bedfile containing location of TSS in column 2')
    parser.add_argument('-radius', metavar='Basepairs', type = int, default = 0,
                        help='Radius to record DS locations. If 0 will take column 2 and 3 as start and end. Otherwise takes column 2 as start, and start + radius as end')
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

    outmat_path = '.'.join([args.outprefix, "mat.gz"])
    outsum_path = '.'.join([args.outprefix, "summary.gz"])

    # Annotate TSS
    bedregions = []
    with gzip.open(args.bedfile, "rt") as inref:
        refreader = csv.reader(inref, delimiter = "\t")
        for row in refreader:
            jchromo = row[0]
            jstrand = row[3]
            jname = row[4]
            if args.radius == 0:
                jstart = int(row[1])
                jend = int(row[2])
            else:
                midpt = int(row[1])
                jstart = midpt - args.radius
                jend = midpt + args.radius + 1
            if jstart < 0:
                # skip ones out of range
                continue
            region = (jchromo, jstart, jend, jstrand, jname)
            bedregions.append(region)


    # Calculate cut positions
    readsthrown = 0
    readstotal = 0
    with pysam.AlignmentFile(args.infile, "rb") as inbam, \
            gzip.open(outmat_path, "wt") as outmat, gzip.open(outsum_path, "wt") as outsum:
        outmatwriter = csv.writer(outmat, delimiter = "\t")
        outsumwriter = csv.writer(outsum, delimiter = "\t")

        for jregion in bedregions:
            jchromo, jstart, jend = jregion[0], jregion[1], jregion[2]
            l = jend - jstart  # length ofvector
            countvec = [0] * l  # init
            # print("Length of countvec:", len(countvec))
            for read in inbam.fetch(contig = jchromo, start = jstart, stop = jend):
                readstotal += 1
                # create output vector
                try:
                    ds = read.get_tag("DS")
                except KeyError:
                    # ds doesnt exist skip 
                    readsthrown += 1
                    continue
                indx = ds - jstart
                if indx < 0 or indx >= len(countvec):
                    # use indx >= len(countvec) because indx == len(countvec) is one number too far past the vector
                    readsthrown += 1
                    continue
                # print(jchromo, jstart, jend)
                # print(indx, ds, jstart)
                # input('debug...')
                # assert indx >= 0 and indx <= len(countvec)
                # print(indx)
                countvec[indx] += 1

            # write row to output file and sumary
            outmatwriter.writerow(countvec)
            region_str = ':'.join([jchromo, '-'.join([str(jstart), str(jend)])])  # chr:123-3456
            outsumwriter.writerow(jregion)
            # print("Wrote reegion", region_str)
            # print(','.join(countvec))

        outmat.close()
        outsum.close()
    print("Reads total:", readstotal)
    print("Reads thrown:", readsthrown)


if __name__ == '__main__':
    main()

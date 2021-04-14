#!/usr/bin/env python
'''
DESCRIPTION

    Calculate dinuc frequency around the cut

FOR HELP

    python calc_dinuc_freq_cuts.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-08-05
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import os
import sys, argparse, datetime
import pysam
from pysamiterators import CachedFasta
from pysam import FastaFile
from singlecellmultiomics.utils import reverse_complement
import numpy as np
import csv



def main():
    parser = argparse.ArgumentParser(description='Calculate dinuc frequency around the cut')
    parser.add_argument('-infile', metavar='INFILE',
                        help='Input bam file')
    parser.add_argument('-outprefix', metavar='OUTPREFIX',
                        help='Outprefix file of each cut and its dinuc frequency (R1 only). Make .txt and .png')
    parser.add_argument('-refpath', metavar='REFFILE',
                        help='Refernce file')
    parser.add_argument('-baseseq', metavar='AGCT sequence', default = 'AT',
                        help='Count sequence as 1 if in baseseq, 0 otherwise')
    parser.add_argument('-mapqthres', metavar='MAPQ', default = 0, type = int, 
                        help='MAPQ default 0')
    parser.add_argument('-upstrm_extend', metavar='distance', default = 0, type = int, 
                        help='Extended left from cut, default 0')
    parser.add_argument('-downstrm_extend', metavar='distance', default = 200, type = int, 
                        help='Extended left from cut, default 200')
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

    inbam = pysam.AlignmentFile(args.infile, "rb")
    reference = CachedFasta(FastaFile(args.refpath))

    readstotal = 0
    readsthrown = 0
    readsbadmapq = 0
    readsbadcigar = 0
    readskept = 0

    dinucfreq_global = [0] * (args.downstrm_extend + args.upstrm_extend)
    dinucfreq_global = np.array(dinucfreq_global)
    print('init')
    print(dinucfreq_global)

    outcsv = args.outprefix + ".csv"
    outsum = args.outprefix + ".summary"
    outfig = args.outprefix + ".png"

    with open(outcsv, "w") as outfile:
        outwriter = csv.writer(outfile, delimiter = "\t")

        for read in inbam.fetch():
            readstotal += 1

            # keep read1 onlyl
            if not read.is_read1:
                readsthrown += 1
                continue

            if read.mapping_quality < args.mapqthres:
                # print("Bad mapping quality, skipping")
                readsthrown += 1
                readsbadmapq += 1
                continue

            chromo = read.reference_name
            chromo_nochr = chromo.replace("chr", '')
            ds = read.get_tag("DS")

            # print(chromo, ds)
            # get ereference 
            if not read.is_reverse:
                strand = "+"
                seq = reference.fetch(chromo_nochr, ds - args.upstrm_extend, ds + args.downstrm_extend)
                # seq = reference.fetch(chromo_nochr, ds, ds + 1)
                # print(seq)
                # print(len(seq))
                # input("testing...")
            else:
                # print("Is reverse:", read.is_reverse)
                # print(chromo, ds, args.downstrm_extend, args.upstrm_extend)
                # print(seq)

                strand = "-"
                seq = reference.fetch(chromo_nochr, ds - args.downstrm_extend, ds + args.upstrm_extend)
                seq = reverse_complement(seq)

                # print(seq)
                # input("Debugging")

            # print("Is reverse:", read.is_reverse)
            # print(seq)
            seq = seq.upper()
            readskept += 1

            # # track dinuc freq
            # dinucfreq = []
            # for s in seq:
            #     dinuc = int(s in args.baseseq)  # 'AT'
            #     dinucfreq.append(dinuc)
            # dinucfreq = np.array(dinucfreq)
            # dinucfreq_global = dinucfreq_global + dinucfreq

            coord = ':'.join([chromo, str(ds), strand])
            outrow = [coord, seq]
            outwriter.writerow(outrow)
    cmd = ' '.join(["gzip", outcsv])
    print(cmd)
    os.system(cmd)

    # write input to plot
    # numpy.savetxt(outsummary, dinucfreq_global, delimiter=",")

    # # make plots
    # xvec = np.array(range(args.upstrm_extend + args.downstrm_extend))
    # plt.plot(xvec, dinucfreq_global / float(readskept), c='r', label = 'output')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(outfig)
    # plt.close()


if __name__ == '__main__':
    main()

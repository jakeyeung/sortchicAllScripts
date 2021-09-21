#!/usr/bin/env python
'''
DESCRIPTION

    Convert refseq chromosome names to ensembl
    From https://www.biostars.org/p/51643/ 
    assembly_report downlaoded fdrom for example: 
    https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_assembly_report.txt

FOR HELP

    python convert_refseq_chromo_to_ucsc.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-08-06
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import gzip, csv

def main():
    parser = argparse.ArgumentParser(description='Convert refseq chromosome names to ensembl')
    parser.add_argument('-infile', metavar='INFILE',
                        help='Input bed-like fille, first column is chromosome name from REFseq eg NC_000067.6')
    parser.add_argument('-assembly_report', metavar='REFFILE',
                        help='Input bed-like fille, first column is chromosome name from REFseq eg NC_000067.6')
    parser.add_argument('-outfile', metavar='OUTFILE',
                        help='Output file, like infile, but first column replaced by new cromosome names')
    parser.add_argument('--add_chr',  action = "store_true",
                        help='Add chromosome prefix toi output')
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

    # Load annots
    refdic = {}
    with open(args.assembly_report, "r") as ref:
        refreader = csv.reader(ref, delimiter = "\t")
        for row in refreader:
            if row[0].startswith("#"):
                continue
            chr_refseq = row[6]  # NC_000067.6 
            chr_ensembl = row[0]  # 1
            if args.add_chr:
                chr_ensembl = ''.join(["chr", chr_ensembl])
            assert chr_refseq not in refdic
            refdic[chr_refseq] = chr_ensembl
    print(refdic)

    with gzip.open(args.infile, "rt") as inf, gzip.open(args.outfile, "wt") as outf:
        jreader = csv.reader(inf, delimiter = "\t")
        jwriter = csv.writer(outf, delimiter = "\t")
        for row in jreader:
            chromo_old = row[0]
            chromo_new = refdic[chromo_old]
            outrow = row
            outrow[0] = chromo_new
            jwriter.writerow(outrow)


if __name__ == '__main__':
    main()

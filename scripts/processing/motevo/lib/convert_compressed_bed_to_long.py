#!/usr/bin/env python
'''
DESCRIPTION

    Convert closest genes (peak to multigene) to long bed (each gene, one line) to make R life easier.

    4th column HMGA1,2.p2;mm10_chr1:3235253-3235753 -> HMGA1,2.p2
    9th column: gene1@gene2 -> each gene separate row
    10th column: dist1@dist2 -> each dist separate row

FOR HELP

    python convert_compressed_bed_to_long.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2016-04-04
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, csv

class BedRow():
    def __init__(self, row, namei = 3, genei = 8, disti = 9, has_dists = True):
        '''
        BedRow with multigene
        '''
        self.row = row
        self.namei, self.genei, self.disti = namei, genei, disti
        self.motifpeak, self.generaw, self.distraw = \
            row[namei], row[genei], row[disti]
        # get motif
        self.motif = self.motifpeak.split(";")[0]
        # convert to lists
        self.genes = self.generaw.split("@")
        if (has_dists):
            if self.distraw == "NA":
                self.dists = ["NA"]
            else:
                self.dists = self.distraw.split("@")  # has commas? need to check
                try:
                    self.dists = [float(i) for i in self.dists]
                except ValueError:
                    # dist has comma eg: 1000,1500. Take min distance
                    for i, dist in enumerate(self.dists):
                        dists_split = [float(d) for d in dist.split(',')]
                        # take min
                        self.dists[i] = min(dists_split)
                # make int
                self.dists = [int(i) for i in self.dists]

    def make_outrow(self, gene, dist = "", motif_only=False):
        '''
        Make outrow, replace gene and dist
        # WARNING: self.row will be modified! original row will be gone
        motif_only: use motif rather than motifpeak
        '''
        outrow = self.row
        if motif_only:
            outrow[self.namei] = self.motif
        else:
            outrow[self.namei] = self.motifpeak
        outrow[self.genei] = gene
        outrow[self.disti] = dist
        return(outrow)

def main():
    parser = argparse.ArgumentParser(description='Convert closest genes (peak to multigene) to long bed (each gene, one line) to make R life easier')
    parser.add_argument('infile', metavar='INFILE',
                        help='multigene bed (from overwrite_nearest_gene_bed.py')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='multigene bed long format suitable for R')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--has_dist', '-d', action='store_true',
                        help='Infile has distance column that needs to be processed')
    parser.add_argument('--motif_only', '-m', action='store_true',
                        help='Take motif only. Default keeps motif and peak together')
    parser.add_argument('--gene_col_i', '-i', metavar='INTEGER', type=int, default=8,
                        help='Index where gene columnn is (python index). Distance_i is 1+i')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    with open(args.infile, 'rb') as readfile, open(args.outfile, 'wb') as writefile:
        bedreader = csv.reader(readfile, delimiter = '\t')
        outwriter = csv.writer(writefile, delimiter = '\t')
        for row in bedreader:
            Row = BedRow(row, genei=args.gene_col_i, disti=args.gene_col_i+1, has_dists = args.has_dist)
            if args.has_dist:
                for gene, dist in zip(Row.genes, Row.dists):
                    outrow = Row.make_outrow(gene, dist)
                    outwriter.writerow(outrow)
            else:
                for gene in Row.genes:
                    outrow = Row.make_outrow(gene)  # makes dist blank
                    outwriter.writerow(outrow, args.motif_only)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
DESCRIPTION

    Parse MARA promoter annotations into a bedfile: NOTE INFILE likely in mm9 need to check

FOR HELP

    python define_promoter_regions.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2016-07-22
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse
import csv

class BedRow():
    def __init__(self, row):
        '''
        Read bed rows, process into an object
        '''
        self.chromo, self.start, self.end = row[0], row[1], row[2]
        self.gene = row[3]
        self.pos = ''.join([self.chromo, ":", self.start, "-", self.end])
        self.posgene = ';'.join([self.gene, self.pos])

    def ProcessPos(self, pos):
        '''chr1:4797368-4798511(+) -> chr1, start, end'''
        pnosuffix = pos[:-3]  # remove (+) or (-)
        chromo = pnosuffix.split(":")[0]
        startend = pnosuffix.split(":")[1]
        start = int(startend.split("-")[0])
        end = int(startend.split("-")[1])
        if start > end:
            # should not happen, complain
            sys.exit("End greater than start. startend: %s" % startend)
        return chromo, str(start), str(end), pnosuffix

    def concat_pos_genename(self, pos, genename, s=";"):
        '''Concatenate pos with genename: chr1:4797368-4798511;Rdh10'''
        return(s.join([pos, genename]))

def main():
    parser = argparse.ArgumentParser(description='Parse MARA promoter annotations into a bedfile')
    parser.add_argument('infile', metavar='INFILE',
                        help='mara_promoters_gene_name_association.bed file')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Bed file of promoter locations. Gene name and position is concatenated.')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
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
        rdr = csv.reader(readfile, delimiter = '\t')
        wtr = csv.writer(writefile, delimiter = '\t')
        for row in rdr:
            Row = BedRow(row)
            wtr.writerow([Row.chromo, Row.start, Row.end, Row.posgene])

if __name__ == '__main__':
    main()

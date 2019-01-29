#!/usr/bin/env python
'''
DESCRIPTION

    Take combined.sites files (output from /Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo/run_merge_motevo_output.sh)
    and convert it to a bed file using the motevo ID.

    Note: the region is a little bit meaningless because I removed "Ns" so I cannot easily be sure whether
    the region is correct. I would need to take the fasta file (with the Ns) and do some scripting
    to extract the exact region, but I'm lazy

    EDIT:
    Being lazy bit me in the butt, but thankfully the Ns do not occur often in the files, so we will record the region as if the software
    ran correctly, then filter out "problematic" IDs downstream

FOR HELP

    python merged_sites_to_bed.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2015-06-03
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse
import os
import csv

class MotRow:
    '''
    Class for handling motevo output row
    '''
    def __init__(self, row):
        self.region = row[0]
        self.strand = row[1]
        self.prob = float(row[2])
        self.motif = row[3]
        self.motevo_id = row[4]

    def merge_row(self, Row):
        '''
        Update parameters by merging with another Row object
        '''
        self.region = ','.join([self.region, Row.region])
        self.strand = ','.join([self.strand, Row.strand])
        self.prob += Row.prob

    def write_to_file(self, writer):
        writer.writerow([self.region, self.strand, self.prob, self.motif, self.motevo_id])

    def get_coord(self):
        '''
        Get bed coordinates from mm10_chr1:84278394-84278894
        '''
        jsplit = self.motevo_id.split("_")[1]
        chromo = jsplit.split(":")[0]
        start = jsplit.split(":")[1].split('-')[0]
        end = jsplit.split(":")[1].split('-')[1]
        return chromo, start, end

    def get_exact_coord(self):
        '''
        Get exact coordinate based on its start and end
        '''
        chromo, start, end = self.get_coord()  # start and end of region of scanned region
        region_start = self.region.split('-')[0]
        region_end = self.region.split('-')[1]
        abs_start = int(start) + int(region_start)
        abs_end = int(start) + int(region_end)
        return chromo, abs_start, abs_end

def get_sites_files(indir):
    dir_files = os.listdir(indir)
    sites_files = []
    for f in dir_files:
        # quick filter for combined.sites
        if f.endswith("combined.sites"):
            fpath = os.path.join(indir, f)
            sites_files.append(fpath)
    return sites_files

def convert_to_bed(inpath, outpath, get_exact = False, add_strand = False):
    with open(inpath, 'rb') as readpath, open(outpath, 'wb') as writepath:
        reader = csv.reader(readpath, delimiter = ' ')
        writer = csv.writer(writepath, delimiter = '\t')
        for row in reader:
            Row = MotRow(row)
            if get_exact:
                chromo, start, end = Row.get_exact_coord()
            else:
                chromo, start, end = Row.get_coord()
            bedid = ';'.join([Row.motif, Row.motevo_id])
            score = Row.prob
            strand = Row.strand
            if not add_strand:
                bedrow = [chromo, start, end, bedid, score]
            else:
                bedrow = [chromo, start, end, bedid, score, strand]
            writer.writerow(bedrow)

def main():
    parser = argparse.ArgumentParser(description='Take combined.sites files (output from /Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo/run_merge_motevo_output.sh) ')
    parser.add_argument('indir', metavar='DIR',
                        help='Input directory containing combined.sites files')
    parser.add_argument('outdir', metavar='DIR',
                        help='Output directory containing combined.sites.bed files')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--get_exact_region', '-e', action='store_true',
                        help='Get exact region of motif rather than just scanned region')
    parser.add_argument('--add_strand', '-s', action='store_true',
                        help='Add strand Felix wants CTCF direction 2017-06-30')
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

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    sites_files = get_sites_files(args.indir)
    for inpath in sites_files:
        infname = os.path.basename(inpath)
        outfname = '.'.join([infname, 'bed'])
        outpath = os.path.join(args.outdir, outfname)
        convert_to_bed(inpath, outpath, get_exact = args.get_exact_region, add_strand = args.add_strand)


if __name__ == '__main__':
    main()

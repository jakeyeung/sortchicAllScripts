#!/usr/bin/env python
'''
DESCRIPTION

    Make count table from bam

FOR HELP

    python make_count_from_bam.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2021-08-26
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import csv


def main():
    parser = argparse.ArgumentParser(description='Make count table from bam')
    parser.add_argument('-inbam', metavar='BAM',
                        help='Input bam')
    parser.add_argument('-metadata_cells', metavar='TXT',
                        help='List of cells to keep')
    parser.add_argument('-metadata_coords', metavar='BED',
                        help='Genome coordinates to look at')
    parser.add_argument('-outfile', metavar='TSV',
                        help='Output count matrix')
    parser.add_argument('-mapq', metavar='MAPQ', type=int, default=40,
                        help='MAPQ cutoff')
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

    # Load cells
    print("Loading cells...")
    cells = []
    with open(args.metadata_cells, "rt") as infmeta:
        metareader = csv.reader(infmeta, delimiter = "\t")
        for row in metareader:
            cells.append(row[0])
    print("Loading cells... done")
    print("Ncells:", len(cells))

    cells_set = set(cells)

    # load genome
    print("Loading genome...")
    coords = []
    coordnames = []
    with open(args.metadata_coords, "rt") as infcoord:
        coordreader = csv.reader(infcoord, delimiter = "\t")
        for row in coordreader:
            # coordtup = (str(row[0]), int(row[1]), int(row[2]))
            startend = "-".join([row[1], row[2]])
            coordtup = ':'.join([row[0], startend])
            coords.append(coordtup)
            if len(row) > 3:
                coordnames.append(row[3])
            else:
                if not args.quiet:
                    print("Warning: no fourth column, appending coord", coordtup)
                coordnames.append(coordtup)

    print("Loading genome... done")
    print("N regions:", len(coords))

    print("Loading cell counts...")
    # iterate bam file, track counts, and location
    countdic = {}

    with pysam.AlignmentFile(args.inbam, "rb") as inbam:
        for coord in coords:
            if coord not in countdic:
                countdic[coord] = {}
                # init dic with zeros for each cell
                for cell in cells:
                    countdic[coord][cell] = 0

            for read in inbam.fetch(region = coord):
                # query name example: TCTCGCGCGCATCGTATGATTGCGGCCATATAGCCT:887823800#3211
                if read.mapping_quality < args.mapq:
                    # skip
                    # print("Skipping because quality too low")
                    # input()
                    continue
                read_barcode = read.query_name.split(":")[0]
                if read_barcode in cells_set:
                    countdic[coord][read_barcode] += 1
    print("Loading cell counts... done")

    print("Writing cell counts...")
    # write count table
    with open(args.outfile, "wt") as outfile:
        outwriter = csv.writer(outfile, delimiter = "\t")
        # write header
        # outheader = ["chromo", "start", "end"]
        outheader = ["coordname"]
        for cell in cells:
            outheader.append(cell)
        outwriter.writerow(outheader)
        for coord, cname in zip(coords, coordnames):
            coordpretty = ';'.join([coord, cname])
            outrow = [coordpretty]
            for cell in cells:
                outrow.append(countdic[coord][cell])
            outwriter.writerow(outrow)
    print("Writing cell counts...")


if __name__ == '__main__':
    main()

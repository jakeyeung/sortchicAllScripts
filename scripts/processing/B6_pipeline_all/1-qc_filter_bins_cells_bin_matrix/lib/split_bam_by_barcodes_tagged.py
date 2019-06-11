#!/usr/bin/env python
'''
DESCRIPTION

    Split bam by barcode specified in a filename

FOR HELP

    python split_bam_by_barcodes.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2018-12-15
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import csv
# from filter_umi_mapq import get_umibc_longheader, get_from_header
import os

def reformat_header(header, add_prefix="chr"):
    new_header = header.to_dict()
    contigs = new_header['SQ']
    new_contigs = []
    for contig in contigs:
        contig['SN'] = ''.join([add_prefix, contig['SN']])
        new_contigs.append(contig)
    new_header['SQ'] = new_contigs
    return pysam.AlignmentHeader.from_dict(new_header)

def scan_bam_get_barcodes(infile, tag="SM"):
    '''
    Read bam file and get barcodes
    '''
    tags = set()
    with pysam.AlignmentFile(infile, "rb") as inbam:
        for read in inbam:
            tags.add(readTag(read, tag))
    return(tags)

def readTag(read, tag, missing='Missing', defective='Defective'):
    '''From /hpc/hub_oudenaarden/bdebarbanson/internalTools/modularDemultiplexer/taggedBamFileToCountTable.py'''
    value=None
    if tag=='chrom':
        return str(read.reference_name)
    if not read.has_tag(tag):
        return missing
    else:
        try:
            value = read.get_tag(tag)
        except Exception as e:
            value = defective
    return value

def main():
    parser = argparse.ArgumentParser(description='Split bam by barcode ')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input bam file')
    parser.add_argument('outdir', metavar='OUTDIR',
                        help='Output directory of bams, split by cell barcodes')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    parser.add_argument('--add_chr_prefix', '-p', action='store_true', help="Add chr prefix to chromosome name")
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

    bc_tag_id = "SM"

    # initialize outbam objects
    bname = os.path.splitext(os.path.basename(args.infile))[0]
    writebamdic = {}  # bam write objs
    sorteddic={}  # output files for sorted bam files
    unsorteddic={}  # tmp files unsorted, will be deleted afterwards
    
    # get barcodes by reading the bam file
    barcodes = scan_bam_get_barcodes(args.infile, tag = bc_tag_id)  # barcodes are sample names
    print(str(len(barcodes)) + " barcodes found in bam file")
    bc_ignored = set()
    with pysam.AlignmentFile(args.infile, "rb") as inbam:
        # change headers?
        if args.add_chr_prefix:
            new_header = reformat_header(inbam.header, add_prefix = "chr")
        for bc in barcodes:
            tmppath = os.path.join(args.outdir, '.'.join([bc, "unsorted", "bam"]))
            outpath = os.path.join(args.outdir, '.'.join([bc, "sorted", "bam"]))
            unsorteddic[bc] = tmppath
            sorteddic[bc] = outpath
            if not args.add_chr_prefix:
                writebamdic[bc] = pysam.AlignmentFile(tmppath, "wb", template = inbam)
            else:
                writebamdic[bc] = pysam.AlignmentFile(tmppath, "wb", header = new_header)
        for read in inbam:
            readbc = readTag(read, bc_tag_id)
            # readbc = get_from_header(read, args.bc_indx)
            if readbc in writebamdic:
                # write to output
                writebamdic[readbc].write(read)
        # close bam files, sort, index and cleanup
        for bc in barcodes:
            writebamdic[bc].close()
            pysam.sort('-o', sorteddic[bc], unsorteddic[bc])
            pysam.index(sorteddic[bc])
        # clean up 
        for bc in barcodes:
            os.remove(unsorteddic[bc])

if __name__ == '__main__':
    main()

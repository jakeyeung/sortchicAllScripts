#!/usr/bin/env python
'''
DESCRIPTION

    Bams need to be filtered by UMI and MAPQ before we do anything...

FOR HELP

    python filter_umi_mapq.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (http://github.com/jakeyeung)
CREATED ON:  2018-12-14
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pysam
import os

def get_umibc(read):
    '''
    ACATAAGAAGA:NS500413:510:H3VGVBGX9:4:23605:24663:8973 -> ACATAAGAAGA
    First 8 nt is UMI, last 3 is cell barcode
    '''
    read_qname = read.query_name
    return(read_qname.split(":")[0])

def get_from_header(read, indx):
    '''
    Is:NS500414;RN:518;Fc:H2GV2BGX9;La:1;Ti:11101;CX:20509;CY:1067;Fi:N;CN:0;aa:CATTTT;aA:CATTTT;aI:35;LY:PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14;RX:CGC;RQ:GGG;BI:185;bc:GTCTGATT;BC:GTCTGATT;QT:GGKKKKKK;MX:NLAIII38 
    -->> GTCTGATT for indx = 16
    -->> CGC for indx = 13
    '''
    jsplit = read.query_name.split(";")[indx]
    return(jsplit.split(":")[1])


def get_umibc_longheader(read, umiindx = 16, bcindx = 13):
    '''
    Long header from Buys?
    Is:NS500414;RN:518;Fc:H2GV2BGX9;La:1;Ti:11101;CX:20509;CY:1067;Fi:N;CN:0;aa:CATTTT;aA:CATTTT;aI:35;LY:PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14;RX:CGC;RQ:GGG;BI:185;bc:GTCTGATT;BC:GTCTGATT;QT:GGKKKKKK;MX:NLAIII38
    4C8U3
    '''
    umi = get_from_header(read, indx = umiindx)
    bc = get_from_header(read, indx = bcindx)
    # concatenate the two
    return(''.join([umi, bc]))

def get_end(read):
    '''
    Return 'R1' if read has positive template_length, R2 otherwise
    '''
    if read.template_length > 0:
        return('R1')
    else:
        return('R2')

def coord_to_bin(coord, window = 1000):
    '''
    6598162 -> neraest 6598000 if window = 1000
    '''
    return(round(coord / window) * window)

def main():
    parser = argparse.ArgumentParser(description='Bams need to be filtered by UMI and MAPQ before we do anything...')
    parser.add_argument('infile', metavar='INFILE',
                        help='bam of a single cell')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='bam filtered by UMI and MAPQ')
    parser.add_argument('outfilesorted', metavar='OUTFILE_SORTED',
                        help='bam filtered by UMI and MAPQ, sorted and indexed')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    parser.add_argument('--mapq_thres', '-m', type=int, default=30,
                        help='Minimum mapping quality. Default 30. bwa output is 0-60?')
    parser.add_argument('--umi_window', '-w', type=int, default=1000,
                        help='When checking UMI, throw out read if things have been read within umi_window basepairs away. Default 1000 to handle paired end?')
    parser.add_argument('--umi_indx', '-i', type=int, default=16,
                        help='umi indx to get from header. Assumes header format of semi-colon separated followed by colon separated of 2 terms')
    parser.add_argument('--bc_indx', '-I', type=int, default=13,
                        help='bc indx to get from header')
    parser.add_argument('--dumpfile', '-d', metavar='DUMPBAM', default = None,
                        help='bam output of UMI duplicates')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    chromos = [''.join(['chr', str(i + 1)]) for i in range(20)] + ['chrX', 'chrY', 'chrM']
    chromos_set = set(chromos)
    bad_chromos = set()
    print(chromos)
    umi_dic = {}  # UMIs are 
    umi_dic_pos = {}  # UMIs positions
    umi_dic_bin = {}  # UMIs by bins, by chromosome
    badreads = 0
    badreadsumi = 0
    goodcounts = 0
    badchromo = 0
    with pysam.AlignmentFile(args.infile, "rb") as bamf:
        with pysam.AlignmentFile(args.outfile, "wb", template=bamf) as outbam:
            if args.dumpfile is not None:
                dumpbam = pysam.AlignmentFile(args.dumpfile, "wb", template=bamf)
            for totalcounts, read in enumerate(bamf):
                if read.mapping_quality < args.mapq_thres:
                    # throw out bad reads
                    badreads += 1
                    continue
                # get UMI-Barcode
                # umibc = get_umibc(read)
                umibc = get_umibc_longheader(read)
                chromo = read.reference_name
                if chromo not in chromos_set:
                    # throw out bad chromo
                    bad_chromos.add(chromo)
                    badchromo += 1
                    continue
                pos = read.reference_start  # left most pos, 0-based
                coord = ':'.join([chromo, str(pos)])
                # get bin within 1kb
                bin = coord_to_bin(pos)
                end = get_end(read)  # Positive or Negative depends on fragment from paired end
                if umibc not in umi_dic:
                    # initialize R1 and R2
                    umi_dic[umibc] = {'R1': 0, 'R2' : 0}  # keep track of reads, they are also indexes!
                    # umi_dic_pos[umibc] = {'R1': [], 'R2' : []}  # keep track of positions
                    umi_dic_bin[umibc] = {'R1': {}, 'R2': {}}  # track bins for UMI counting
                    for c in chromos:
                        for end in ['R1', 'R2']:
                            umi_dic_bin[umibc][end][c] = set()
                else:
                    # update umi_dics
                    # if umibc counted in same bin in same end (R1 or R2) then it's bad
                    if bin in umi_dic_bin[umibc][get_end(read)][chromo]:
                        # umi_dic_pos[umibc][get_end(read)].append(coord)  # record the bad read 
                        umi_dic[umibc][get_end(read)] += 1  # record duplicate reads 
                        # already in bin, then don't add it
                        badreadsumi += 1
                        dumpbam.write(read)
                        continue
                # reads here are unique (within a window) and high quality, write them to outbam and move on
                # umi_dic_pos[umibc][get_end(read)].append(coord)  # only append if it's a bad read 
                umi_dic_bin[umibc][get_end(read)][chromo].add(bin)
                goodcounts += 1
                outbam.write(read)
            if args.dumpfile is not None:
                dumpbam.close()
    # sort and index bam
    pysam.sort('-o', args.outfilesorted, args.outfile)
    pysam.index(args.outfilesorted)
    # remove temporarily file
    os.remove(args.outfile)

    # Print arguments supplied by user. Ideally as log file because of the \n 
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c')); print('\n')
        print('Command line inputs:'); print('\n')
        print(CMD_INPUTS); print('\n')
        print ('Argparse variables:'); print('\n')
        print(ARG_INPUTS); print('\n')
        # list total counts
        print("Reads from these chromosomes thrown out: %s" % bad_chromos)
        print("Total counts: %s" % totalcounts); print('\n')
        print("High quality, unique counts: %s" % goodcounts); print('\n')
        print("Bad quality counts: %s" % badreads); print('\n')
        # write duplicate UMI counts as a table
        for umibc in umi_dic:
            if umi_dic[umibc]['R1'] > 1 and umi_dic[umibc]['R2'] > 1:
                # print('%s\t%s\t%s\t%s\%s\n' %\
                #         (umibc, umi_dic[umibc]['R1'], 
                #             umi_dic[umibc]['R2'], 
                #             ','.join(umi_dic_pos[umibc]['R1']),
                #             ','.join(umi_dic_pos[umibc]['R2'])))
                print('%s\t%s\t%s\n' %\
                        (umibc, umi_dic[umibc]['R1'], 
                            umi_dic[umibc]['R2']))

if __name__ == '__main__':
    main()

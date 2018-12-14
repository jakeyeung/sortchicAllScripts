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

def get_umibc(read):
    '''
    ACATAAGAAGA:NS500413:510:H3VGVBGX9:4:23605:24663:8973 -> ACATAAGAAGA
    First 8 nt is UMI, last 3 is cell barcode
    '''
    read_qname = read.query_name
    return(read_qname.split(":")[0])

def get_end(read):
    '''
    Return 'R1' if read has positive template_length, R2 otherwise
    '''
    if read.template_length > 0:
        return('R1')
    else:
        return('R2')

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
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    umi_dic = {}  # UMIs are 
    umi_dic_pos = {}  # UMIs are 
    badreads = 0
    badreadsumi = 0
    goodcounts = 0
    with pysam.AlignmentFile(args.infile, "rb") as bamf:
        with pysam.AlignmentFile(args.outfile, "wb", template=bamf) as outbam:
            for totalcounts, read in enumerate(bamf):
                if read.mapping_quality < args.mapq_thres:
                    # throw out bad reads
                    badreads += 1
                    continue
                # get UMI-Barcode
                umibc = get_umibc(read)
                chromo = read.reference_name
                pos = read.reference_start  # left most pos, 0-based
                coord = ':'.join([chromo, str(pos)])
                end = get_end(read)  # Positive or Negative depends on fragment from paired end
                if umibc not in umi_dic:
                    # initialize R1 and R2
                    umi_dic[umibc] = {'R1': 0, 'R2' : 0}  # keep track of reads, they are also indexes!
                    umi_dic_pos[umibc] = {'R1': [], 'R2' : []}  # keep track of positions
                else:
                    # update umi_dics
                    # if umibc already in dic, can be counted if R1 or R2 are zero, 
                    # otherwise it's a badreadsumi
                    if umi_dic[umibc][get_end(read)] > 1:  # > not >= b/c we added 1 alrdy
                        umi_dic[umibc][get_end(read)] += 1
                        umi_dic_pos[umibc][get_end(read)].append(coord)
                        badreadsumi += 1
                        continue
                # reads here are unique and high quality, keep write them to outbam and move on
                umi_dic[umibc][get_end(read)] += 1
                umi_dic_pos[umibc][get_end(read)].append(coord)
                goodcounts += 1
                outbam.write(read)
    # sort and index bam
    pysam.sort('-o', args.outfilesorted, args.outfile)
    pysam.index(args.outfilesorted)

    # Print arguments supplied by user
    if not args.quiet:
        if args.logfile is not None:
            sys.stdout = open(args.logfile, "w+")
        print(datetime.datetime.now().strftime('Code output on %c')); print('\n')
        print('Command line inputs:'); print('\n')
        print(CMD_INPUTS); print('\n')
        print ('Argparse variables:'); print('\n')
        print(ARG_INPUTS); print('\n')
        # list total counts
        print("Total counts: %s" % totalcounts); print('\n')
        print("High quality, unique counts: %s" % goodcounts); print('\n')
        print("Bad quality counts: %s" % badreads); print('\n')
        # write duplicate UMI counts as a table
        for umibc in umi_dic:
            if umi_dic[umibc]['R1'] > 1 and umi_dic[umibc]['R2'] > 1:
                print('%s\t%s\t%s\t%s\%s\n' %\
                        (umibc, umi_dic[umibc]['R1'], 
                            umi_dic[umibc]['R2'], 
                            ','.join(umi_dic_pos[umibc]['R1']),
                            ','.join(umi_dic_pos[umibc]['R2'])))

if __name__ == '__main__':
    main()

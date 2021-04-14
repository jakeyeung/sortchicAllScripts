import sys
import argparse
from Bio import SeqIO
import re
import gzip

def main():
    parser = argparse.ArgumentParser(description='Parse gff')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input gff')
    # parser.add_argument('outfile', metavar='OUTFILE',
    #                     help='Output txt')
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

    # # Print arguments supplied by user
    # if not args.quiet:
    #     if args.logfile is not None:
    #         sys.stdout = open(args.logfile, "w+")
    #     print(datetime.datetime.now().strftime('Code output on %c'))
    #     print('Command line inputs:')
    #     print(CMD_INPUTS)
    #     print ('Argparse variables:')
    #     print(ARG_INPUTS)
        
    count = 0
    filename = args.infile

    for seq_record in SeqIO.parse(filename, 'gb'):
        # print(seq_record)
        pattern = re.compile("^NC_(\d){6}.(\d)+")
        if (pattern.match(seq_record.id)):
            for feature in seq_record.features:
                # print(feature)
                if (feature.type == 'mRNA'):
                    strd = feature.location.strand
                    start = feature.location.nofuzzy_start
                    end = feature.location.nofuzzy_end
                    if (strd == 1):
                        strand = '+'
                    else:
                        strand = '-'
                    gene = feature.qualifiers['gene'][0]
                    description = feature.qualifiers['product'][0]
                    transcript_id = feature.qualifiers['transcript_id'][0]
                    db_xref = feature.qualifiers['db_xref']
                    print('%s\t%s\t%s\t%s\tmRNA\t%s\t%s\t%s' % (seq_record.id, start, end, strand, gene, transcript_id, description))
                if (feature.type == 'CDS'):
                    strd = feature.location.strand
                    start = feature.location.nofuzzy_start
                    end = feature.location.nofuzzy_end
                    if (strd == 1):
                        strand = '+'
                    else:
                        strand = '-'
                    gene = feature.qualifiers['gene'][0]
                    description = 'n/a'
                    if 'protein_id' in feature.qualifiers:
                        transcript_id = feature.qualifiers['protein_id'][0]
                    else:
                        transcript_id = 'n/a'
                    db_xref = 'n/a'
                    print('%s\t%s\t%s\t%s\tCDS\t%s\t%s\t%s' % (seq_record.id, start, end, strand, gene, transcript_id, description))
                if (feature.type == 'exon'):
                    strd = feature.location.strand
                    start = feature.location.nofuzzy_start
                    end = feature.location.nofuzzy_end
                    if (strd == 1):
                        strand = '+'
                    else:
                        strand = '-'
                    gene = feature.qualifiers['gene'][0]
                    transcript_id = "n/a"
                    description = "n/a"
                    db_xref = 'n/a'
                    print('%s\t%s\t%s\t%s\texon\t%s\t%s\t%s' % (seq_record.id, start, end, strand, gene, transcript_id, description))


if __name__ == '__main__':
    main()

# filename = "/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.gff"
# count = 0



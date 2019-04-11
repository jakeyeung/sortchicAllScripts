#!/usr/bin/env python
'''
DESCRIPTION

    Convert npz numpy zip object to textfile

FOR HELP

    python npz_to_textfile.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2019-03-28
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime, os
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Convert npz numpy zip object to textfile')
    parser.add_argument('infile', metavar='INFILE',
                        help='npz file')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='output tab delimited file')
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

    assert os.path.exists(args.infile)
    bname, _ = os.path.splitext(args.infile)
    dat = np.load(args.infile)
    jheader = '\t'.join(dat['labels'])
    dat_mat = dat['matrix']
    with open(args.outfile, 'wb') as f:
        np.savetxt(f, dat_mat, delimiter = "\t", header = jheader, comments = "")

if __name__ == '__main__':
    main()

#!/usr/bin/env python
'''
DESCRIPTION

    Coupling matrix to file numpy ndarray to file

FOR HELP

    python downstream_scot.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-11-24
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import os
import sys, argparse, datetime
import numpy as np
import pickle

def main():
    parser = argparse.ArgumentParser(description='Coupling matrix to file')
    parser.add_argument('infile', metavar='INFILE',
                        help='Input transport coupling matrix pickle')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='Output .txt file of transport coupling matrix')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    parser.add_argument('--logfile', '-l', metavar='LOGFILE', default = None,
                        help='Write arguments to logfile')
    args = parser.parse_args()

    # # store command line arguments for reproducibility
    # CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # # store argparse inputs for reproducibility / debugging purposes
    # args_dic = vars(args)
    # # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]  # for python2
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]  # for python3
    # ARG_INPUTS = ' '.join(ARG_INPUTS)

    # # Print arguments supplied by user
    # if not args.quiet:
    #     if args.logfile is not None:
    #         sys.stdout = open(args.logfile, "w+")
    #     print(datetime.datetime.now().strftime('Code output on %c'))
    #     print('Command line inputs:')
    #     print(CMD_INPUTS)
    #     print ('Argparse variables:')
    #     print(ARG_INPUTS)

    with open(args.infile, 'rb') as infobj:
        gamma = pickle.load(infobj)

    print("Writing to file")
    np.savetxt(args.outfile, gamma, delimiter = "\t")
    # gamma.tofile(args.outfile, sep = "\t")
    # print("Compressing file")
    # os.system("gzip %s" %args.outfile)

if __name__ == '__main__':
    main()

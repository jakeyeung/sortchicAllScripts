#!/usr/bin/env python
'''
DESCRIPTION

    Run SCOT

FOR HELP

    python run_scot_on_scchic_script.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2020-11-23
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import pickle
import sys
# sys.path.insert(0, "/usr/local/lib/python3.7/site-packages")
# sys.path.insert(0, "/Users/yeung/projects/unbalanced_gromov_wasserstein")
sys.path.insert(0, "/hpc/hub_oudenaarden/jyeung/code_for_analysis/unbalanced_gromov_wasserstein")
sys.path.insert(0, "/hpc/hub_oudenaarden/jyeung/code_for_analysis/SCOT")
import ot
import src.utils as ut
import src.evals as evals
from src.scot import *
# from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
# from sklearn.decomposition import PCA
import numpy as np 
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Run SCOT')
    parser.add_argument('-infile1', metavar='INFILE',
                        help='Cell by latent features matrix')
    parser.add_argument('-infile2', metavar='INFILE',
                        help='Cell by latent features matrix, second matrix')
    parser.add_argument('-outprefix', metavar='OUTFILE',
                        help='Pickle output prefix. Will add coupling and log to the output')
    parser.add_argument('-knn', metavar='Integer', type=int, default=50,
                        help='K nearest neighbors')
    parser.add_argument('-eps', metavar='Float', type=float, default=0.0005,
                        help='Epsilon for entropy penalty')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    args = parser.parse_args()

    # # store command line arguments for reproducibility
    # CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # # store argparse inputs for reproducibility / debugging purposes
    # args_dic = vars(args)
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]
    # ARG_INPUTS = ' '.join(ARG_INPUTS)

    # # Print arguments supplied by user
    # if not args.quiet:
    #     print(datetime.datetime.now().strftime('Code output on %c'))
    #     print('Command line inputs:')
    #     print(CMD_INPUTS)
    #     print ('Argparse variables:')
    #     print(ARG_INPUTS)

    inf1 = args.infile1
    inf2 = args.infile2
    outf1 = ''.join([args.outprefix, "_coupling.pkl"])
    outf2 = ''.join([args.outprefix, "_logoutput.pkl"])
    print(outf1)
    print(outf2)
    X = pd.read_csv(inf1, sep = '\t')
    Y = pd.read_csv(inf2, sep = '\t')
    k = args.knn
    e = args.eps

    gamma, logoutput = scot(X, Y, k, e, mode="connectivity", metric="correlation", returnCoupling = True, balanced = False)

    pickle.dump(gamma, outf1)
    pickle.dump(logoutput, outf2)
    

if __name__ == '__main__':
    main()

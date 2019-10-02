#!/usr/bin/env python
'''
DESCRIPTION

    Run LDA 

FOR HELP

    python run_LDA_gensim.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2019-07-10
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import time
import sys, argparse, datetime
import gensim
from gensim import corpora
from gensim.test.utils import datapath
import scipy.io
import logging
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
import numpy as np
import csv

def main():
    parser = argparse.ArgumentParser(description='Run LDA ')
    parser.add_argument('infile', metavar='INFILE',
                        help='mm input')
    parser.add_argument('outfile', metavar='OUTFILE',
                        help='LDA output')
    parser.add_argument('ntopics', metavar='INTEGER', type=int, help='Number of topics')
    # parser.add_argument('--eta', metavar='ETA', type=float, help='Eta prior', default = 0.1)
    parser.add_argument('--eta', metavar='ETA', help='Eta prior', default = 0.1)
    parser.add_argument('--alpha', metavar='ALPHA', help='Alpha prior', default = 'bytopics')
    parser.add_argument('--etashift', metavar='ETA SHIFT', type=float, help='Add constant to eta', default = 0)  # for variational bayes
    parser.add_argument('--alphashift', metavar='ALPHA SHIFT', type=float, help='Add constant to alpha', default = 0)  # for variational bayes
    parser.add_argument('--ncores', metavar='NCORES', type=int, help='Number of cores', default = 1)
    parser.add_argument('--npasses', metavar='PASSES', type=int, help='Number of passes through the data', default = 100)
    parser.add_argument('--chunksize', metavar='PASSES', type=int, help='Chunk size', default = 384)  # chunk approximately plates
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

    scipy_sparse_matrix=scipy.io.mmread(args.infile)
    corpus = gensim.matutils.Sparse2Corpus(scipy_sparse_matrix)
    # corpus = corpora.MmCorpus(args.infile)
    # run LDA
    if args.alpha == "bytopics":
        jalpha = 50./args.ntopics
        alphavec = args.ntopics * [jalpha]
    elif args.alpha == "auto":
        alphavec = "auto"  # https://radimrehurek.com/gensim/models/ldamodel.html
    else:
        # assume jalpha is a float
        jalpha = float(args.alpha)
        alphavec = args.ntopics * [jalpha]
    # shift priors for variational bayes: https://arxiv.org/pdf/1205.2662.pdf
    if isinstance(alphavec, list):
        alphavec = [x + args.alphashift for x in alphavec]
    if args.eta == "auto":
        jeta = "auto"
    else:
        jeta = float(args.eta)
        # shift
        jeta = jeta + args.etashift

    if jeta == "auto":
        print("Using autotuning for eta")
    if alphavec == "auto":
        print("Using autotuning for alpha")

    if args.ncores > 1:
        print("Running LdaMulticore...")
        lda = gensim.models.ldamulticore.LdaMulticore(corpus=corpus, num_topics=args.ntopics, alpha=alphavec, eta=jeta, minimum_probability = 0, chunksize=args.chunksize, passes=args.npasses, workers=args.ncores)
    else:
        print("Running LdaModel single core...")
        lda = gensim.models.ldamulticore.LdaModel(corpus=corpus, num_topics=args.ntopics, update_every=0, alpha=alphavec, eta=jeta, minimum_probability=0, chunksize=args.chunksize, passes=args.npasses)

    lda.save(args.outfile)

    all_topics = lda.get_document_topics(corpus, per_word_topics=False, minimum_probability=0, minimum_phi_value=0)
    termtopmat = lda.get_topics()

    termtop_outf = ''.join([args.outfile, ".topic_to_terms.csv"])
    doctop_outf = ''.join([args.outfile, ".doc_to_topics.csv"])
    # write topic to terms
    np.savetxt(termtop_outf, termtopmat, delimiter="\t")    

    # write doc to topics
    with open(doctop_outf, mode="w") as wf:
        writer = csv.writer(wf, delimiter = "\t")
        for d in all_topics:
            # each d is a list of tuples of (topic_id, topic_weight)
            outrow = [i[1] for i in d]
            assert(abs(sum(outrow) - 1) < 0.00001)  # check things sum up to ~1
            writer.writerow(outrow)

if __name__ == '__main__':
    main()

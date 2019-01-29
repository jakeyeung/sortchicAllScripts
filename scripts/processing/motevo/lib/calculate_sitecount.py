#! /software/bin/python

import csv
import sys
import os
import argparse
import pysam
import numpy as np
import re
# from Bio import SeqIO

PSEUDO_COUNT = 1.0
INDEX = {'A':0, 'C':1, 'G':2, 'T':3}
INDEX_INV = {0:'A', 1:'C', 2:'G', 3:'T'}
LOG025 = np.log(0.25)


def arguments():
    parser = argparse.ArgumentParser(description="""
    Given a set of WMs and a FASTA file, this code calculates the
    sitecount matrix for the sequences.
    """)
    parser.add_argument('-w', '--wm', dest='wmDir', action='store',
                       type=str, required=True,
                       help='The input directory that contains WM files')
    parser.add_argument('-f', '--fasta', dest='fastafile', action='store',
                       type=str, required=True,
                       help='The input FASTA file.')
    parser.add_argument('-c', '--cutoff', dest='cutoffDir', action='store',
                       type=str, required=True,
                       help='The directory which holds the cutoff values.')
    parser.add_argument('-C', '--CutoffValue', dest='jcutoff', type = float, default = 0.1)
    parser.add_argument('-E', '--EMprior', dest='EMprior', type = int, default = 0)  # 0 or 1
    parser.add_argument('-B', '--BGprior', dest='BGprior', type = float, default = 0.9995)  # 0 to 1
    parser.add_argument('-o', '--output', dest='outputDir', action='store',
                       type=str, required=True,
                       help='The name of the output directory.')
    args = parser.parse_args()
    return args


def loadWM(wmFile):
    """ returns the the normalized weight matrix in LOG space """
    content = [l for l in open(wmFile) if re.search('^\d+\s+\d+\.*\d*\s+\d+\.*\d*', l)]
    wm = np.zeros(len(content)*4).reshape(len(content),4)
    for c, i in zip(content, range(len(content))):
        wm[i] = map(lambda x: float(x)+PSEUDO_COUNT, c.split()[1:5])
    row_sums = np.sum(wm, axis=1)
    return np.log(wm / row_sums[:, np.newaxis])


def giveRevWM(wm):
    return wm[::-1, ::-1]


def calculateMaxScore(wm):
    score = np.sum([np.max(row) for row in wm])
    return score - wm.shape[0]*LOG025


def informationContent(wm):
    info = 0.
    for row in wm:
        for elem in row:
            info -= np.exp(elem)*(elem - LOG025)
    return info


def getWMLength(wm):
    return wm.shape[0]


def wmScore(seq, wm):
    return np.sum([wm[i][INDEX[s]] for s, i in zip(seq, xrange(wm.shape[0]))]) \
      - wm.shape[0]*LOG025


def calculateSitecount(seq, wm, wm_r, maxscore):
    length = wm.shape[0]
    motifcount = 0.
    for i in xrange(len(seq) - length):
        site = seq[i:i+length]
        if not re.search('N+', site):
            motifcount += np.exp(wmScore(site, wm) - maxscore)
            motifcount += np.exp(wmScore(site, wm_r) - maxscore)
    return motifcount


def scanInputSequences(fastafile, WMs, WMs_rev, WMMaxScores, outf):
    with open(fastafile) as inf:
        for record in SeqIO.parse(inf, "fasta"):
            seqid = str(record.id)
            seq = str(record.seq).upper()
            motifsCounts = []
            for i in xrange(len(WMs)):
                motifsCounts.append(calculateSitecount(seq, WMs[i], \
                                                       WMs_rev[i], WMMaxScores[i]))
            outf.write('\t'.join([
                    seqid,
                    '\t'.join(map(lambda x: '%.6f' % x, motifsCounts)) + '\n',
                    ]))


def clusterJob(path, wmfile, wmname, fastafile, outputDir, cutoffDir, param_file):
    fname = os.path.join(path, wmname + '.sh')
    with open(fname, "w") as outf:
        outf.write('\n'.join([
            "#!/bin/bash",
            "#$ -S /bin/bash",
            "#$ -M j.yeung@hubrecht.eu",
            "#$ -m beas",
            "#$ -N %s" % wmname,
            "#$ -l h_vmem=8000000",
            "#$ -l h_rt=2:00:00",
            "",
            ]))
        cmd = ' '.join([
            '/hpc/hub_oudenaarden/jyeung/software/motevo_ver1.11/bin/motevo',
            fastafile,
            '"%s"' % param_file,
            '"%s"\n' % wmfile,
            ])
        outf.write(cmd)
    # os.system('bsub < "%s"' % fname)
    # os.system('bsub < "%s"' % fname)
    # os.system("qsub %s" % fname)  # dont run it yet
    return fname


def main():
    args = arguments()
    # path = '/scratch/cluster/monthly/somidi/Tissue.Specificity/data/sitecounts'
    path = ''
    WMNames = sorted(os.listdir(args.wmDir))
    WMFiles = [os.path.join(args.wmDir, wm) for wm in WMNames]
    for wmfile, wmname in zip(WMFiles, WMNames):
        with open("%s.param" % wmname , "w") as param_file:
            param_file.write("\n".join([
                "refspecies mm10",
                "TREE (mm10:1);",
                "",
                "Mode TFBS",  # different modes, use TFBS
                # "EMprior 0",  # Use a constant prior (do not update it)
                "EMprior %s" % (args.EMprior),
                "",
                "markovorderBG 0",
                # "bgprior 0.9995",
                "bgprior %s" % (args.BGprior),
                "bg A 0.25",
                "bg T 0.25",
                "bg C 0.25",
                "bg G 0.25",
                "",
                "restrictparses 0",
                "priordiff 0.001",
                "sitefile %s.sites" % (os.path.join(args.outputDir, wmname)),
                "priorfile %s.prior" % (os.path.join(args.outputDir,  wmname)),
                "minposterior %s\n" % args.jcutoff,  # put filter greatly reduces output size
                ]))

        clusterJob(path, wmfile, wmname, args.fastafile, args.outputDir, args.cutoffDir, "%s.param" % wmname)
        # exit()


if __name__ == '__main__':
    main()

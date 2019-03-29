#!/bin/sh
# Jake Yeung
# bam_to_bigwig.sh
# coverage bam to bigwig 
# 2019-03-20

# effective genomesize https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# tutorial https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#usage-example-for-chip-seq

inbam=$1
outbw=$2

[[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
[[ -e $outbw ]] && echo "$outbw found, exiting for safety" && exit 1
bamCoverage --bam $inbam -o $outbw \
    --binSize 10
    --normalizeUsing RPGC
    --effectiveGenomeSize 2652783500

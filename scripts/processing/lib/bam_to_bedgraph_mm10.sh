#!/bin/sh
# Jake Yeung
# bam_to_bedgraphsh
# coverage bam to bedgraph, so we can do logtransform afterwards
# 2019-03-20

# effective genomesize https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# tutorial https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html#usage-example-for-chip-seq

inbam=$1
outbedgraph=$2

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

[[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
[[ -e $outbedgraph ]] && echo "$outbedgraph found, exiting for safety" && exit 1
bamCoverage --bam $inbam -o $outbedgraph --binSize 100 --normalizeUsing CPM -of bedgraph

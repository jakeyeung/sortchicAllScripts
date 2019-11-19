#!/bin/sh
# Jake Yeung
# 0-index_bwa_zebrafish.sh
#  
# 2019-10-31

inf="/hpc/hub_oudenaarden/group_references/ensembl/98/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

bwa index $inf

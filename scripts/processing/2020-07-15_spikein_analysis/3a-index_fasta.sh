#!/bin/sh
# Jake Yeung
# 3a-index_fasta.sh
#  
# 2020-07-15

jmem='64G'
jtime='6:00:00'

infasta="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92_WithLambdaPhage.fa"

BNAME="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/nohupsout/BwaIndex.out"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3; bwa index $infasta"

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=BwaIndex --wrap "$cmd"

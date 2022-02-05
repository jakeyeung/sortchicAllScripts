#!/bin/sh
# Jake Yeung
# 8-run_cellranger.sh
#  
# 2021-06-16

jmem='32G'
jtime='2:00:00'

bs="/hpc/hub_oudenaarden/jyeung/software/cellranger-atac-2.0.0/bin/cellranger-atac"
ref="/hpc/hub_oudenaarden/jyeung/data/databases/CellRanger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz"

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/fastqs_check"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/rename"

# cd $indir

# f1=${indir}/SRR12638101_1.fastq.gz
# f2=${indir}/SRR12638101_2.fastq.gz

id="SRR12638101"
cmd="cd $indir; cellranger-atac count --id ${id} --reference $ref --fastqs ${indir}"

BNAME=${indir}/${id}.sbatch_log
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${id} --wrap "$cmd"

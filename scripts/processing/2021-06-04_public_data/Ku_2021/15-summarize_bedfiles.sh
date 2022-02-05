#!/bin/sh
# Jake Yeung
# 15-summarize_bedfiles.sh
#  
# 2021-06-15

jmem='4G'
jtime='1:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output/summarized"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/Ku_2021/summarize_bedfiles.R"

for f in `ls -d $indir/*.bed.gz`; do
    fbase=$(basename $f)
    fbase=${fbase%%.*}
    outf=${outdir}/${fbase}.summarized_by_cell_barcodes.bed

    BNAME=${outdir}/${fbase}.sbatch.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $f $outf"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $f $outf
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${fbase} --wrap "$cmd"
done

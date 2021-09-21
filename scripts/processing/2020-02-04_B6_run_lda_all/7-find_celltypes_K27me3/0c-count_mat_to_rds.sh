#!/bin/sh
# Jake Yeung
# 0c-count_mat_to_rds.sh
#  
# 2020-08-29

jmem='8G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/convert_txt_to_rds_count_table.R"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix/countTables.unfixed"
outdir=${indir}

for f in `ls -d $indir/*.csv.gz`; do
    fbase=$(basename $f)
    fbase=${fbase%.*}
    fbase=${fbase%.*}

    BNAME=${outdir}/${fbase}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${fbase}.rds
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $f $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${fbase} --wrap "$cmd"
done


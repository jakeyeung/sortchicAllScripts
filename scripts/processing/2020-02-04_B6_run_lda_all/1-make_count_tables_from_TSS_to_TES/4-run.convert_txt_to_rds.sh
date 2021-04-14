#!/bin/sh
# Jake Yeung
# 4-run.convert_txt_to_rds.sh
#  
# 2020-08-23

jmem='8G'
jtime='2:00:00'

# WRAP UP
while [[ `squeue -u jyeung | grep count | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/convert_txt_to_rds_count_table.R"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40/countTablesAndRZr1only.NewFilters"
outdir="${indir}/rds_objs"
[[ ! -d $outdir ]] && mkdir $outdir

for inf in `ls -d $indir/*50kb.csv`; do
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatch
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outfile=${outdir}/${bname}.rds
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outfile -format bed"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=Convert_${bname} --wrap "$cmd"
done


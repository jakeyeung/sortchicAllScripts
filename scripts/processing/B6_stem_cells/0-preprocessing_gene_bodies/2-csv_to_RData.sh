#!/bin/sh
# Jake Yeung
# 2-csv_to_rds.sh
#  
# 2019-06-24

jmem='8G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic_gastru/scripts/processing/lib/csv_to_RData_TSS.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/count_tables_genebodies"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3700/countTables_genebodies"
outdir=$indir
for inf in `ls -d $indir/*.csv`; do
    # inf="/hpc/hub_oudenaarden/jyeung/data/scchic_gastru/histone-mods/count_tables_merged_TSS/${prefix}.filtered.mapq_60.TSS_50000.v2pipeline.csv"
    # [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    # outf="/hpc/hub_oudenaarden/jyeung/data/gastru_scchic/histone-mods/count_tables_merged_TSS"
    # outf="/hpc/hub_oudenaarden/jyeung/data/scchic_gastru/histone-mods/count_tables_merged_TSS/mm_H3K36me3.filtered.mapq_60.TSS_50000.withname.bugfixed.rds"
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=$outdir/$bname.csvToRData
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=$outdir/$bname.RData
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
done
# wait

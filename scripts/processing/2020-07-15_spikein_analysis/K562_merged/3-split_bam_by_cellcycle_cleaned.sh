#!/bin/sh
# Jake Yeung
# 3-split_bam_by_cellcycle_cleaned.sh
#  
# 2020-10-18

jmem='16G'
jtime='12:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
mapq=40

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K562_clean_by_cellcycle"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_${mapq}"
outdir="$indir/bams_CellCycleSorted_split_by_cellcycle"
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1  # make dir before hand

for jmark in $jmarks; do
    bname="K562_CellCycleSorted_${jmark}.merged.sorted.tagged"
    inf=${indir}/${bname}.bam
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    # annotname="BM_AllMerged.${jmark}.cell_cluster_table.txt"
    annotname="K562_cleaned_by_cellcycle.${jmark}.txt"
    annotfile=${annotdir}/${annotname}
    [[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix"
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N split_${jmark} -m beas -M j.yeung@hubrecht.eu
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

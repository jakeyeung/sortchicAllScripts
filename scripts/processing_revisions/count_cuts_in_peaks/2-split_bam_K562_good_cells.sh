#!/bin/sh
# Jake Yeung
# split_bam_K562_good_cells.sh
#  
# 2022-01-10

# WRAP UP

while [[ `squeue -u jyeung | grep merge | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='12:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
jmarks="k4me1 k4me3 k27me3 k9me3"
mapq=40

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/K562/merged_bams"
annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/metadata/K562"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1

outdir="$indir/bams_split_by_good_cells_NoChrPrefix"
[[ ! -d $outdir ]] && mkdir $outdir

for jmark in $jmarks; do
    # indir="${inmain}/K562_${jmark}"
    bname="K562_revisions_no_spikeins.${jmark}"
    inf=${indir}/${bname}.bam
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    # annotname="K562_cleaned_by_G1filt.${jmark}.txt"
    annotname="metadata_for_splitting.K562_${jmark}.0.8_0.5_3000.2022-01-10.txt"
    annotfile=${annotdir}/${annotname}
    [[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq}"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

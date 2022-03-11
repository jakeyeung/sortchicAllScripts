#!/bin/sh
# Jake Yeung
# 1-merge_bams.sh
# Merge bam before splitting
# 2022-01-10

jmem='16G'
jtime='24:00:00'


inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/K562"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/K562/merged_bams"

jmarks="k4me1 k4me3 k27me3 k9me3"
cell="K562"

for jmark in $jmarks; do
    indir="${inmain}/${cell}_${jmark}"
    inbams="${indir}/*.bam"
    echo $inbams
    outbam=${outdir}/${cell}_revisions_no_spikeins.${jmark}.bam
    BNAME=${outdir}/${cell}_${jmark}_sbatch_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate SCMO_2021; samtools merge $outbam $inbams; samtools index $outbam"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=merge_${jmark} --wrap "$cmd"
done

#!/bin/sh
# Jake Yeung
# 4b-merge_bams_by_marks.sh
#  
# 2020-08-11

jmem='32G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562/tagged_bams/merged_bams"

# merge H3K4me1
# jmark="H3K4me1"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do

    BNAME=${outdir}/${jmark}.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    infs="${indir}/K562-EtOH-${jmark}-G1*.sorted.tagged.bam"
    outf="${outdir}/K562-EtOH-${jmark}.merged.sorted.tagged.bam"

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $infs"
    echo $cmd

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=merge_${jmark} --wrap "$cmd"
done


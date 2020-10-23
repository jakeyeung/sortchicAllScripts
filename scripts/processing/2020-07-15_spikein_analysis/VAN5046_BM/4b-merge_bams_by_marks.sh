#!/bin/sh
# Jake Yeung
# 4b-merge_bams_by_marks.sh
#  
# 2020-09-12

jmem='16G'
jtime='3:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir=${indir}/merged_bams
[[ ! -d $outdir ]] && mkdir $outdir

jmarks="H3K4me1 H3K4me3 H3K27me3"



cd $indir
for jmark in $jmarks; do
    BNAME=${outdir}/${jmark}.sbatch.output
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outf="${outdir}/PZ-ChIC-mouse-BM-${jmark}-merged.sorted.tagged.bam"
    infs=`ls -d ${indir}/PZ-ChIC-mouse-BM-${jmark}*.sorted.tagged.bam | tr '\n' ' '`
    echo $jmark
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $infs"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
done


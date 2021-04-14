#!/bin/sh
# Jake Yeung
# 4-sort_index_tag_bam.sh
#  
# 2019-12-19

jmem='16G'
jtime='6:00:00'

inf1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams/PZ-ChIC-mouse-BM-H3K27me3-merged.sorted.tagged.bam"
inf2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5230_BM/tagged_bams/merged_bams/PZ-ChIC-BM-rep3-H3K27me3-merged.sorted.tagged.bam"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merged_across_runs"
outname="PZ-ChIC_H3K27me3_merged.VAN5046_VAN5230"
outf="${outdir}/${outname}.sorted.tagged.bam"

BNAME=${outdir}/${outname}.log
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $inf1 $inf2; samtools index $outf"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MergeBams --wrap "$cmd"



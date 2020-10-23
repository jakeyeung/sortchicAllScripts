#!/bin/sh
# Jake Yeung
# 3-merge_bams_by_mark.sh
#  
# 2020-09-24

jmem='32G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams"
outdir="${indir}/merged_by_mark"

jmarks="k36me3 k27me3 k9me3"

for jmark in $jmarks; do

    BNAME=${outdir}/${jmark}_merge.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

   outf="${outdir}/CG-ChIC-TAPS-RPE-exp1-${jmark}.tagged.bam"
   infs="$indir/*${jmark}*.bam"
   # echo $infs
   cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $infs"
   # echo "samtools merge $outf $indir/*${jmark}.bam"
   sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
done

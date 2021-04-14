#!/bin/sh
# Jake Yeung
# 4b-merge_bams.sh
#  
# 2020-10-07

jmem='32G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams"
outdir="${indir}/merged_bams"

# cd $indir

# merge K4me1

jmark="H3K4me3"
infs=$(ls -d $indir/*rep3-${jmark}*.bam | tr '\n' ' ')
outf=${outdir}/PZ-BM-rep3-${jmark}.merged.sorted.tagged.bam
# echo $infs

[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

BNAME=$outdir/${jmark}_merge.log
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf"
echo $cmd
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${outf}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_merge --wrap "$cmd"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams"

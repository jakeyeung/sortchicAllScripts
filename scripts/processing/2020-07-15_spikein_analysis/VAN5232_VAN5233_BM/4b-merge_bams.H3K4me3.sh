#!/bin/sh
# Jake Yeung
# 4b-merge_bams.sh
#  
# 2020-10-07

jmem='32G'
jtime='6:00:00'

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM/mouse/tagged_bams"
outdir="${indir}/merged_bams"
[[ ! -d $outdir ]] && mkdir $outdir

# cd $indir

# merge K4me1
jmark="H4K4me1"  # typo
jmarkout="H3K4me1"
infs=$(ls -d $indir/PZ-BM-rep2-${jmark}*.bam | tr '\n' ' ')
outf=${outdir}/PZ-BM-rep2-${jmarkout}.merged.sorted.tagged.bam
# echo $infs

BNAME=$outdir/H3K4me1_merge.log
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf"
echo $cmd
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${outf}_%j.log --ntasks=1 --nodes=1 --job-name=H3K4me1_merge --wrap "$cmd"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM/mouse/tagged_bams"

#!/bin/sh
# Jake Yeung
# 7-merge_tagged_bams_by_plate_H3K4me1.sh
# 
# 2019-12-14

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks"
[[ ! -d $outdir ]] && mkdir $outdir

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

# Merge H3K4me3: Unenriched, Linneg, and StemCells
bin=$indir/*H3K4me3*.retagged.bam
bout=$outdir/H3K4me3-BM_Linneg_SC-merged.tagged.bam
samtools merge -f $bout $bin --output-fmt BAM
samtools index $bout


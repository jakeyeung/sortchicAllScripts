#!/bin/sh
# Jake Yeung
# 8a-merge_tagged_bams_by_plate.sh
# Create one tagged bam for each mark 
# 2019-11-26

indir=/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/merge_by_mark"
[[ ! -d $outdir ]] && mkdir $outdir

# jmarks="H3K4me1 H3K4me3 H3K9me3"
jmarks="H3K9me3"
# jmark="H3K4me1"


for jmark in $jmarks; do
    infsBM=$indir/PZ-ChIC-ZFWKM-${jmark}*.retagged.bam
    infsSC=$indir/PZ-ChIC-ZFWKMCD41plus-${jmark}-*.retagged.bam
    outbam=$outdir/${jmark}-WKM_CD41-merged.tagged.retagged.bam
    [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    samtools merge -f $outbam $infsBM $infsSC --output-fmt BAM
    samtools index $outbam
done


#!/bin/sh
# Jake Yeung
# 8a-merge_tagged_bams_by_plate.sh
# Create one tagged bam for each mark 
# 2019-11-26
# H3K27me3 has no ZF enrichment yet

indir=/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/merge_by_mark"
[[ ! -d $outdir ]] && mkdir $outdir

# jmarks="H3K4me1 H3K4me3 H3K9me3"
jmarks="H3K27me3"
# jmark="H3K4me1"


for jmark in $jmarks; do
    infsBM=$indir/PZ-ChIC-ZFWKM-${jmark}*.retagged.bam
    outbam=$outdir/${jmark}-WKM-merged.tagged.retagged.bam
    [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    samtools merge -f $outbam $infsBM --output-fmt BAM
    samtools index $outbam
done


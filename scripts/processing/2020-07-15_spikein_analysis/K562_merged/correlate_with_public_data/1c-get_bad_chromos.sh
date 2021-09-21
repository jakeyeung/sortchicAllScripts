#!/bin/sh
# Jake Yeung
# 1c-get_bad_chromos.sh
#  
# 2020-11-22

inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt/K562_AllMerged_H3K4me1.merged.sorted.tagged.G1filt.sorted.bam"

outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/badchromos.txt"

samtools view -H $inbam | awk '(NR == 24) || (NR >= 27 && NR <= 288) { print $2 }' | cut -d":" -f2 > $outf

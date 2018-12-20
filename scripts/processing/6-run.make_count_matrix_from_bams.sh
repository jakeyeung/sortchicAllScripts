#!/bin/sh
# Jake Yeung
# 6-run.make_count_matrix_from_bams.sh
# Make count matrix from bams using R 
# 2018-12-19

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $workdir

# rs="/home/hub_oudenaarden/jyeung/projects/scChiC/processing/make_count_matrix_from_bams.R"
rs="processing/make_count_matrix_from_bams.R"
peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.blacklistfilt.broadPeak"

jchip="H3K4me1"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1

# get paths to bams
bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc"

tmpf="/tmp/JY_${jchip}_bamlist.out"

[[ -e $tmpf ]] && echo "$tmpf must not already exist" && exit 1

for b in $(ls -d $bmain/PZ*$jchip*/*.sorted.bam); do
    echo $b >> $tmpf
done

echo "Temp bam files found in $tmpf"

# now run Rscript
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats/PZ-BM-${jchip}.merged.mat"

Rscript $rs $tmpf $peakf $outf

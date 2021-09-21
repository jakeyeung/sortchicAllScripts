#!/bin/sh
# Jake Yeung
# 1b-put_bigwig_one_dir.sh
#  
# 2020-11-22

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_ENCODE_rerun_from_Buys"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bigwigs_G1filt_split_by_G1filt.for_chipseq_comparison"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/bigwigs_to_compare"
[[ ! -d $outdir ]] && mkdir $outdir

for b in `ls -d $indir1/*.bw`; do
    ln -s $b $outdir
done

for b in `ls -d $indir2/*.bw`; do
    ln -s $b $outdir
done

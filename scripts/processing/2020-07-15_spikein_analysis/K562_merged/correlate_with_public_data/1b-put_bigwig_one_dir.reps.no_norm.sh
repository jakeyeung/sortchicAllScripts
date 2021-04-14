#!/bin/sh
# Jake Yeung
# 1b-put_bigwig_one_dir.sh
#  
# 2020-11-22

# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_ENCODE_rerun_from_Buys"
# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/chipseq_from_ENCODE_first_submission"
# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/chipseq_from_ENCODE_first_submission/bsize_1000"
# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562"
# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_ENCODE_first_submission"
inbase="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams"
# indir1="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K4me1"
# indir2="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K4me3"
# indir3="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K27me3"
# indir4="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K9me3"
indir5="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bigwigs_G1filt_split_by_G1filt.for_chipseq_comparison"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/bigwigs_to_compare.reps.no_norm"
# [[ -d $outdir ]] && echo "$outdir found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
for mark in $jmarks; do
    indirtmp=${inbase}/${mark}
    for b in `ls -d $indirtmp/*10kb.bw`; do
        bname=$(basename $b)
        bname=${bname%.*}
        echo $bname
        ln -s $b $outdir/${mark}.${bname}.bw
    done
done

for b in `ls -d $indir5/*.bw`; do
    echo $b
    ln -s $b $outdir
done

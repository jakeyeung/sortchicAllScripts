#!/bin/sh
# Jake Yeung
# compare_multi_bigwigs.Server.sh
# Compare multiple big wigs 
# 2019-03-21

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/renamed"

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
ctypes="Neu Pro"

for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison"
    [[ ! -d $outdir ]] && mkdir $outdir
    for mark in $marks; do
        echo $mark
        infbase="$inbase/${mark}_${ctype}.bw"
        [[ ! -e $infbase ]] && echo "$infbase not found, exiting" && exit 1

        chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_2019-03-20"  # 9 clusters
        infcompares=$(ls -d $chicbwdir/${mark}_cluster_*.bw)

        echo "multiBigwigSummary bins -b $infbase $infcompares -o $outdir/${mark}_${ctype}_comparison.npz --binSize 100000"
        multiBigwigSummary bins -b $infbase $infcompares -o $outdir/${mark}_${ctype}_comparison.npz --binSize 100000 
        ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done



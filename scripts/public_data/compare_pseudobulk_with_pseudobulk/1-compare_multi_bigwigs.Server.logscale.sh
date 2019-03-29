#!/bin/sh
# Jake Yeung
# compare_multi_bigwigs.Server.sh
# Compare multiple big wigs 
# 2019-03-21

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_2019-03-20/log1p_bigwigs"  # 9 clusters
[[ ! -d $chicbwdir ]] && echo "$chicbwdir not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
ctypes="Multimark"

# H3K4me1 vs H3K4me3, H3K27me3, and H3K9me3
markref="H3K4me1"
marks="H3K4me3 H3K27me3 H3K9me3"

infblist=$(ls -d $chicbwdir/${markref}_cluster_*.bw)
infblist=$(echo $infblist | tr '\n' ' ')

suffix="log1p_"

for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_${suffix}comparison"
    [[ ! -d $outdir ]] && mkdir $outdir
    for mark in $marks; do
        echo $mark
        infclist=$(ls -d $chicbwdir/${mark}_cluster_*.bw)
        infclist=$(echo $infclist | tr '\n' ' ')
        outf="$outdir/${mark}_${ctype}_${markref}_vs_${mark}_${suffix}comparison.npz"
        echo "multiBigwigSummary bins -b $infblist $infclist -o $outf --binSize 100000"
        multiBigwigSummary bins -b $infblist $infclist -o $outf --binSize 100000
        ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done


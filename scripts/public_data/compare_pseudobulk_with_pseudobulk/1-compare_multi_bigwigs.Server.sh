#!/bin/sh
# Jake Yeung
# compare_multi_bigwigs.Server.sh
# Compare multiple big wigs 
# 2019-03-21

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_2019-03-20"  # 9 clusters

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
ctypes="Multimark"

# H3K4me1 vs H3K4me3, H3K27me3, and H3K9me3
markref="H3K4me1"
marks="H3K4me3 H3K27me3 H3K9me3"

infblist=$(ls -d $chicbwdir/${markref}_cluster_*.bw)
infblist=$(echo $infblist | tr '\n' ' ')

binsize="1000"

n=0
maxjobs=3

for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison"
    [[ ! -d $outdir ]] && mkdir $outdir
    for mark in $marks; do
        echo $mark
        infclist=$(ls -d $chicbwdir/${mark}_cluster_*.bw)
        infclist=$(echo $infclist | tr '\n' ' ')
        outf="$outdir/${mark}_${ctype}_${markref}_vs_${mark}_comparison_binsize-${binsize}.npz"
        [[ -e $outf ]] && echo "$outf found, continuing" && continue
        echo "multiBigwigSummary bins -b $infblist $infclist -o $outf --binSize $binsize"
        multiBigwigSummary bins -b $infblist $infclist -o $outf --binSize $binsize&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
        # ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done



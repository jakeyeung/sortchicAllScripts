#!/bin/sh
# Jake Yeung
# 1-compare_multi_bigwigs.Server.build95.sh
# Compare bigwigs with build95 
# 2019-03-28


# compare neutrophils H3K4me1 versus all H3K4me1 clusters

suffix="build95_B6"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}"

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
ctypes="Multimark"

# H3K4me1 vs H3K4me3, H3K27me3, and H3K9me3
# markref="H3K4me3"
# markref="H3K4me3"
# marks="H3K4me1 H3K27me3 H3K9me3"

# markref="H3K4me1"
# marks="H3K4me3 H3K27me3 H3K9me3"
# 

markref="H3K27me3"
marks="H3K4me1 H3K4me3 H3K9me3"

infblist=$(ls -d $chicbwdir/${markref}_cluster_*.bw)
infblist=$(echo $infblist | tr '\n' ' ')

binsize="100000"

n=0
maxjobs=6

for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_binsize-${binsize}_comparison_${suffix}"
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
wait

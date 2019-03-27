#!/bin/sh
# Jake Yeung
# compare_multi_bigwigs.Server.sh
# Compare multiple big wigs 
# 2019-03-21

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs/renamed"

marks="H3K4me1 H3K4me3 H3K27me3"
ctypes="MatBcell HSC ProB"

for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison"
    [[ ! -d $outdir ]] && mkdir $outdir
    for mark in $marks; do
        if [[ "$mark" == "H3K27me3" && "$ctype" == "HSC" ]]; then
            continue
        fi
        # echo $mark
        for infbase in `ls -d $inbase/${mark}_${ctype}*.bw`; do
            rep=$(basename $infbase)
            rep=${rep%%.*}
            rep=$(echo $rep | cut -d"_" -f2)
            chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_2019-03-20"  # 9 clusters
            infcompares=$(ls -d $chicbwdir/${mark}_cluster_*.bw)
            infcompares=$(echo $infcompares | tr '\n' ' ')
            # echo $infcompares
            # echo $infcompares
            # echo "multiBigwigSummary bins -b $infbase $infcompares -o $outdir/${mark}_${ctype}_comparison.npz --binSize 100000"
            # multiBigwigSummary bins -b $infbase $infcompares -o $outdir/${mark}_${ctype}_${rep}_comparison.npz --binSize 100000 
            # echo "multiBigwigSummary bins -b $infbase $infcompares -o $outdir/${mark}_${ctype}_${rep}_comparison.npz --binSize 100000"
            outf="$outdir/${mark}_${ctype}_${rep}_comparison.npz" 
            [[ -e $outf ]] && echo "$outf  found, continuing" && continue
            echo "multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000"
            multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000

            ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

        done

    done
done



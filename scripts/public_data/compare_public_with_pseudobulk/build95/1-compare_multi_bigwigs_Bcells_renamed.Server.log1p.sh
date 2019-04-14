#!/bin/sh
# Jake Yeung
# compare_multi_bigwigs.Server.sh
# Compare multiple big wigs 
# 2019-03-21

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs/renamed_relative_to_input"

marks="H3K4me1 H3K4me3 H3K27me3"
ctypes="MatBcell HSC ProB"
suffix="build95"
subdir="log1p_bigwigs"

chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}_2019-03-28/${subdir}"
[[ ! -d $chicbwdir ]] && echo "$chicbwdir not found, exiting" && exit 1

blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
[[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1

n=0
maxjobs=1
for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison_${suffix}_${subdir}"
    [[ ! -d $outdir ]] && mkdir $outdir
    for mark in $marks; do
        if [[ "$mark" == "H3K27me3" && "$ctype" == "HSC" ]]; then
            continue
        fi
        # echo $mark
        for infbase in `ls -d $inbase/${mark}_${ctype}*.bw`; do
            rep=$(basename $infbase)
            rep=${rep%%.*}
            rep=$(echo $rep | cut -d"_" -f3)
            infcompares=$(ls -d $chicbwdir/${mark}_cluster_*.bw)
            infcompares=$(echo $infcompares | tr '\n' ' ')
            outf="$outdir/${mark}_${ctype}_${rep}_comparison.npz"
            [[ -e $outf ]] && echo "$outf  found, continuing" && continue
            echo "multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000 -bl $blist"
            multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000 -bl $blist
            if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            	# define maxjobs and n using maxjobsn skeleton
                wait # wait until all have finished (not optimal, but most times good enough)
                echo $n wait
            fi
            # ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
        done

    done
done
echo "Done script"

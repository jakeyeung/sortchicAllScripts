#!/bin/sh
# Jake Yeung
# 1-compare_multi_bigwigs_Lara-Astiaso.Server.sh
# Compare multiple big wigs 
# 2019-03-26

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

jdate="2019-03-28"
suffix="build95"
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}_${jdate}"  # 9 clusters
[[ ! -d $chicbwdir ]] && echo "$chicbwdir not found, exiting" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/renamed"
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/Lara-Astiaso_2014_${suffix}"

[[ ! -d $outbase ]] && mkdir $outbase

marks="H3K4me1 H3K4me3"  # active marks only 

cd $inbase  # so ls doesn't give full name
ctypes=$(ls -d *.bw | sed 's/.bw//g' | cut -d"_" -f2 | tr '\n' ' ')

blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
[[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1

echo $ctypes

n=0
maxjobs=8
for ctype in $ctypes; do
    echo $ctype
    for mark in $marks; do
        # echo $mark
        infbase=$inbase/${mark}_${ctype}.bw
        [[ ! -e $infbase ]] && echo "$infbase not found, exiting" && exit 1
        infcompares=$(ls -d $chicbwdir/${mark}_cluster_*.bw)
        infcompares=$(echo $infcompares | tr '\n' ' ')
        outf="$outbase/${mark}_${ctype}_comparison.npz" 
        [[ -e $outf ]] && echo "$outf  found, continuing" && continue
        # echo "multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000"
        multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000 -bl $blist&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
        
        # ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done

echo "Done"

#!/bin/sh
# Jake Yeung
# compare_multi_bigwigs.Server.InputNorm.Eryth.sh
# Compare multiple big wigs 
# 2019-03-21

# compare neutrophils H3K4me1 versus all H3K4me1 clusters

suffix="build95"
subdir="log1p_bigwigs"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_GenomeRes_2014/bigwig_mm10/renamed_relative_to_input"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1

chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}_2019-03-28/$subdir"
[[ ! -d $chicbwdir ]] && echo "$chicbwdir not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
ctypes="Erythrobl Megakaryo"  # erythro has H3K9me3, Megakaryo does not
reps="rep1 rep2"

blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
[[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1

for ctype in $ctypes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison_${suffix}_${subdir}"
    [[ ! -d $outdir ]] && mkdir $outdir
    for mark in $marks; do
      if [[ $ctype == "Megakaryo" && $mark == "H3K9me3" ]]; then
        echo "skipping $ctype and $mark"
      	continue
      fi  
        echo $mark
        for rep in $reps; do
            infbase="$inbase/${mark}_${ctype}_${rep}_reltoinput.bw"
            [[ ! -e $infbase ]] && echo "$infbase not found, exiting" && exit 1

            infcompares=$(ls -d $chicbwdir/${mark}_cluster_*.bw)
            infcompares=$(echo $infcompares | tr '\n' ' ')

            # echo $infcompares
            outf="$outdir/${mark}_${ctype}_${rep}_comparison.npz"
            echo "multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000 -bl $blist"
            multiBigwigSummary bins -b $infbase $infcompares -o $outf --binSize 100000 -bl $blist
            ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
        done
    done
done
echo "Done script"


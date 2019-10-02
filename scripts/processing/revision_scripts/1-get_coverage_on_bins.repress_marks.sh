#!/bin/sh
# Jake Yeung
# 1-get_coverage_on_bins.sh
# Get .bam file and make .bigwig coverage from a bedfile
# 2019-08-12



inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/revision_objs/BM_bins.2019-08-12.bed"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

suffix="build95_B6"  # do not do stringent
outsuff="_revisions"
chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}"  # H3K4me3 clusters
[[ ! -d $chicbwdir ]] && echo "$chicbwdir not found, exiting" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
# inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/renamed"
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_GenomeRes_2014/bigwig_mm10/renamed"
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/Wu_GenomeRes_2014_${suffix}${outsuff}"

[[ ! -d $outbase ]] && mkdir $outbase

marks="H3K27me3 H3K9me3"  # do all

# ctypes=$(ls -d *.bw | sed 's/.bw//g' | cut -d"_" -f2 | tr '\n' ' ')

# blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
# [[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1

n=0
maxjobs=4
for mark in $marks; do
    cd $inbase  # so ls doesn't give full name
    ctypes=$(ls -d ${mark}*.bw | sed 's/.bw//g' | cut -d"_" -f2 | sort | uniq | tr '\n' ' ')
    for ctype in $ctypes; do
        echo $ctype
        infbases=$(ls -d $inbase/${mark}_${ctype}_*.bw)
        # echo $infbases
        for infbase in $infbases; do
            bbase=$(basename $infbase)
            bbase=${bbase%.*}
            rep=$(echo ${bbase} | cut -d"_" -f3)
            echo $rep
            infcompares=$(ls -d $chicbwdir/${mark}_cluster_1.bw)  # compare with just first cluster
            infcompares=$(echo $infcompares | tr '\n' ' ')
            outf="$outbase/${mark}_${ctype}_${rep}_comparison.npz"
            outftxt="$outbase/${mark}_${ctype}_${rep}_comparison.tab"
            [[ -e $outf ]] && echo "$outf  found, continuing" && continue
            multiBigwigSummary BED-file -b $infbase $infcompares -o $outf --outRawCounts $outftxt --BED $inf&
            if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                    # define maxjobs and n using maxjobsn skeleton
                wait # wait until all have finished (not optimal, but most times good enough)
                echo $n wait
            fi
        done
    done
done

echo "Done"



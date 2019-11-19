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
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/renamed"
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/Lara-Astiaso_2014_${suffix}${outsuff}"

[[ ! -d $outbase ]] && mkdir $outbase

marks="H3K4me1 H3K4me3"  # do all

cd $inbase  # so ls doesn't give full name
ctypes=$(ls -d *.bw | sed 's/.bw//g' | cut -d"_" -f2 | tr '\n' ' ')

# blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
# [[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1

echo $ctypes
n=0
maxjobs=4
for ctype in $ctypes; do
    echo $ctype
    for mark in $marks; do
        # echo $mark
        infbase=$inbase/${mark}_${ctype}.bw
        [[ ! -e $infbase ]] && echo "$infbase not found, exiting" && exit 1
        infcompares=$(ls -d $chicbwdir/${mark}_cluster_1.bw)  # compare with just first cluster
        infcompares=$(echo $infcompares | tr '\n' ' ')
        outf="$outbase/${mark}_${ctype}_comparison.npz"
        outftxt="$outbase/${mark}_${ctype}_comparison.tab"
        [[ -e $outf ]] && echo "$outf  found, continuing" && continue
        multiBigwigSummary BED-file -b $infbase $infcompares -o $outf --outRawCounts $outftxt --BED $inf&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
        # ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done

echo "Done"



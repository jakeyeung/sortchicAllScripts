#!/bin/sh
# Jake Yeung
# 1-call_hidden_domains_on_clustered_bams.sh 
# 2019-04-15

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

suffix="_from_traj"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_build95_B6${suffix}"
[[ ! -d $outmain ]] && mkdir $outmain

marks="H3K4me1"

inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_build95_B6${suffix}"

n=0
maxjobs=12
for mark in $marks; do
    for inf in $(ls -d $inbase/${mark}_cluster_*.bam); do
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        bname=$(basename $inf)
        bname=${bname%.*}

        minlength=1000  # seems to work well at 1k?
        bname=$bname.$minlength.cutoff
        [[ ! -d $outmain ]] && mkdir $outmain
        outdir=$outmain/$bname

        [[ -d $outdir ]] && echo "$outdir found, continuing" && continue

        [[ ! -d $outmain ]] && mkdir $outmain
        [[ ! -d $outdir ]] && mkdir $outdir
        BNAME=$outdir/$bname.log

        . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done
wait

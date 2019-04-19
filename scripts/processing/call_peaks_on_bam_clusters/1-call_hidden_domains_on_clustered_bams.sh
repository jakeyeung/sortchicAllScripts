#!/bin/sh
# Jake Yeung
# 1-call_hidden_domains_on_clustered_bams.sh 
# 2019-04-15

jmem='4G'
jtime='1:00:00'

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_build95"
[[ ! -d $outmain ]] && mkdir $outmain

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_build95_2019-03-28"

for mark in $marks; do
    for inf in $(ls -d $inbase/${mark}_cluster_*.bam); do
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        bname=$(basename $inf)
        bname=${bname%.*}

        minlength=1000  # seems to work well at 1k?
        bname=$bname.$minlength.cutoff
        [[ ! -d $outmain ]] && mkdir $outmain
        outdir=$outmain/$bname

        # [[ -d $outdir ]] && echo "$outdir found, skipping $mark" && continue

        [[ ! -d $outmain ]] && mkdir $outmain
        [[ ! -d $outdir ]] && mkdir $outdir
        BNAME=$outdir/$bname.log

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
        # echo "hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname"
    done
done


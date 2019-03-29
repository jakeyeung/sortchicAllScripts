#!/bin/sh
# Jake Yeung
# 2c-call_hidden_domains_on_merged_bams.sh
#  
# 2019-01-08

jmem='4G'
jtime='1:00:00'

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output_build95"
[[ ! -d $outmain ]] && mkdir $outmain

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/BM_${mark}_merged.bam"
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/add_chromo_prefix"

# # run on GPU node??
# n=0
# maxjobs=4

for mark in $marks; do
    # echo $mark

    inf=$inbase/BM_${mark}_merged.bam
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}

    minlength=1000  # seems to work well at 1k?
    bname=$bname.$minlength.cutoff
    [[ ! -d $outmain ]] && mkdir $outmain
    outdir=$outmain/$bname

    [[ -d $outdir ]] && echo "$outdir found, skipping $mark" && continue

    [[ ! -d $outmain ]] && mkdir $outmain
    [[ ! -d $outdir ]] && mkdir $outdir
    BNAME=$outdir/$bname.log


    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname&
    # if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    #     # define maxjobs and n using maxjobsn skeleton
    #     wait # wait until all have finished (not optimal, but most times good enough)
    #     echo $n wait
    # fi
done


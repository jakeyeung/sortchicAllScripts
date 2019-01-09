#!/bin/sh
# Jake Yeung
# 2c-call_hidden_domains_on_merged_bams.sh
#  
# 2019-01-08

jmem='16G'
jtime='3:00:00'

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.txt"

[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"

for mark in $marks; do
    # echo $mark
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM_${mark}_merged.bam"
    bname=$(basename $inf)
    bname=${bname%.*}

    minlength=1000  # seems to work well at 1k?
    bname=$bname.$minlength.cutoff
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output"
    [[ ! -d $outmain ]] && mkdir $outmain
    outdir=$outmain/$bname

    [[ -d $outdir ]] && echo "$outdir found, skipping $mark" && continue

    [[ ! -d $outmain ]] && mkdir $outmain
    [[ ! -d $outdir ]] && mkdir $outdir
    BNAME=$outdir/$bname.log

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
done


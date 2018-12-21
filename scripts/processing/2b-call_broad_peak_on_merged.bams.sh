#!/bin/sh
# Jake Yeung
# 2b-call_broad_peak_on_merged.bams.sh
# Call broad peak on merged bams 
# 2018-12-18

jmem='10G'
jtime='3:00:00'

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM_H3K4me1_merged.bam"
bname=$(basename $inf)
bname=${bname%.*}

cutoff=0.3
minlength=1111
bname=$bname.$cutoff.$minlength.cutoff
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output"
outdir=$outmain/$bname

[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outdir ]] && mkdir $outdir
BNAME=$outdir/$bname.log

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2; macs2 callpeak --cutoff-analysis -t $inf --broad -g mm --pvalue $cutoff --broad-cutoff $cutoff --min-length $minlength -n $bname --outdir $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err

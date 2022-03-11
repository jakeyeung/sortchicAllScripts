#!/bin/sh
# Jake Yeung
# 3-run_hiddendomains_on_CUTnTAG.sh
#  
# 2021-06-05

jmem='4G'
jtime='1:00:00'

jmark="H3K4me2"

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.hg38.filt.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K4me2/merged_beds"

outmain="${inmain}/hiddendomains_output"
[[ ! -d $outmain ]] && mkdir $outmain

maxcounts=10
minlength=1000  # seems to work well at 1k?


inf="${inmain}/K562_${jmark}_20181120_allmerged.bed"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

bname=$(basename $inf)
bname=${bname%.*}
bname=$bname.$minlength.cutoff
outdir=$outmain/$bname

[[ ! -d $outdir ]] && mkdir $outdir
BNAME=$outdir/$bname.log

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname -B"

sbatch --time=${jtime} --mem-per-cpu=${jmem} --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${mark}_${bname} --wrap "$cmd"


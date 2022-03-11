#!/bin/sh
# Jake Yeung
# 14-run_hiddendomains_on_bam.sh
# Wu
# 2021-06-21

jmem='4G'
jtime='1:00:00'

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.hg38.filt.nochr.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/merged_beds_dedup"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1_notrim_bowtie2"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outmain="${inmain}/hiddendomains_output"
[[ ! -d $outmain ]] && mkdir $outmain

maxcounts=10
minlength=2500  # seems to work well at 1k?
mapq=0

inf="${inmain}/SRR12638101_1.sorted.bam"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

bname=$(basename $inf)
bname=${bname%.*}
bname=$bname.$minlength.cutoff
outdir=$outmain/$bname

[[ ! -d $outdir ]] && mkdir $outdir
BNAME=$outdir/$bname.log

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname -q $mapq"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"


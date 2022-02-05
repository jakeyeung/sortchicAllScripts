#!/bin/sh
# Jake Yeung
# 14-run_hiddendomains_on_bam.sh
# Wu
# 2021-06-21

jmem='4G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/run_hidden_domains.R"

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.hg38.filt.nochr.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/merged_beds_dedup"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1_notrim_bowtie2"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outmain="${inmain}/hiddendomains_output"
[[ ! -d $outmain ]] && mkdir $outmain

# maxcounts=10
mapq=0
maxcount="25"
mincount="1"
minprob="0.6"
minlength="1000"


inf="${inmain}/SRR12638101_1.sorted.bam"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

bname=$(basename $inf)
bname=${bname%.*}
bname=$bname.$minlength.cutoff
outdir=$outmain/$bname

[[ ! -d $outdir ]] && mkdir $outdir
BNAME=$outdir/$bname.log

dname=${bname}
inf="${outmain}/${dname}/${dname}_treatment_bins.txt"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
outf="${outmain}/${dname}/${dname}_domains.txt"
[[ -e $outf ]] && echo "$outf found, continuing" && continue

chromonames="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

# cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname -q $mapq"
cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -maxreadcount $maxcount -minreadcount $mincount -minprob $minprob -chromonames ${chromonames}"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"


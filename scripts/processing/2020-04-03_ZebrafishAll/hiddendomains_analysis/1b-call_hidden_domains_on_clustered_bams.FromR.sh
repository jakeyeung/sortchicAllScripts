#!/bin/sh
# Jake Yeung
# 1-call_hidden_domains_on_clustered_bams.sh 
# 2019-04-15

jmem='4G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/run_hidden_domains.R"

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/danRer11.chrom.sizes"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs.FromR"
[[ ! -d $outmain ]] && mkdir $outmain
# [[ ! -d $inmain ]] && mkdir $inmain

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"


maxcount="20"
mincount="-10"
minprob="0.6"
chromonames="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25"

for indir in `ls -d $inmain/PZ*`; do
    dname=$(basename $indir)
    outdir=${outmain}/${dname}
    [[ ! -d $outdir ]] && mkdir $outdir

    inf="${indir}/${dname}_treatment_bins.txt"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outf="${outdir}/${dname}_domains.txt"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    BNAME=$outdir/$dname.log

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -maxreadcount $maxcount -minreadcount $mincount -minprob $minprob -chromonames ${chromonames}" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N "R_HD_$dname"
done


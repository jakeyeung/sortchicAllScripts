#!/bin/sh
# Jake Yeung
# 5-call_peaks.sh
# Call peaks pseudobulk 
# 2020-10-29

jmem='4G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/run_hidden_domains.R"

# chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/danRer11.chrom.sizes"
chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.hg38.filt.nochr.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs"
inmain="pwd/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt_NoChrPrefix/hiddendomains_outputs"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
# outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs.FromR"
outmain="pwd/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt_NoChrPrefix/hiddendomains_outputs.FromR"
[[ ! -d $outmain ]] && mkdir $outmain
# [[ ! -d $inmain ]] && mkdir $inmain

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

maxcount="20"
mincount="-10"
minprob="0.6"

chromonames="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

for indir in `ls -d $inmain/K562*`; do
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


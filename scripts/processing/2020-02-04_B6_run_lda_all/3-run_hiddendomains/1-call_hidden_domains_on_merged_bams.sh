#!/bin/sh
# Jake Yeung
# 1-call_hidden_domains_on_merged_bams.sh
# 2020-02-14

jmem='16G'
jtime='3:00:00'

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.filt.txt"
[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
minlength=1000  # seems to work well at 1k?
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output"
[[ ! -d $outbase ]] && mkdir $outbase


for mark in ${marks}; do
    outmain=${outbase}/hd_clusters.${mark}.minlength_${minlength}
    [[ ! -d $outmain ]] && mkdir $outmain
    for inf in `ls -d ${indir}/${mark}*.sorted.bam`; do
        bname=$(basename $inf)
        bname=${bname%.*}
        bname=${bname}.minlength_${minlength}
        outdir=$outmain/$bname
        [[ -d $outdir ]] && echo "$outdir found, skipping $mark" && continue
        [[ ! -d $outdir ]] && mkdir $outdir
        BNAME=$outdir/$bname.qsub
        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname"
        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N ${bname}
    done
done

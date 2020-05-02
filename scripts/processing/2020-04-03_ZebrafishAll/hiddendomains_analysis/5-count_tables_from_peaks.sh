#!/bin/sh
# Jake Yeung
# 3-count_tables_from_peaks.sh
#  
# 2020-04-24


jmem='32G'
jtime='6:00:00'

mapq=40
# dist=
ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.test.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# blacklist="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.stringent.2020-04-14.corrfilt/correlated_bins.bed"
blacklist="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-14.corrfilt/correlated_bins.chromofixed.bed"

jsuffix="imputevarfilt.lessstringent.mapq_40"
bamdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jsuffix}.remerged_by_marks"
[[ ! -d $bamdir ]] && echo "$bamdir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.HiddenDomains.${jsuffix}"
[[ ! -d $outdir ]] && mkdir $outdir

for jmark in $jmarks; do
    inbed="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs.FromR/hd_merged.${jmark}.minlength_1000/merged.${jmark}.minlength_1000.cutoff_analysis.merged.withchr.annotated.bed"
    [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

    bname="${jmark}.${jsuffix}"
    inbamname="${bname}.remerged.bam"
    inbamenamei="${bname}.remerged.bam.bai"

    inbam=${bamdir}/${inbamname}
    outf=${outdir}/${bname}.countTable.HiddenDomains.csv

    BNAME=${outdir}/bamToCountTable_qsub_${bname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed -blacklist $blacklist"
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed -blacklist $blacklist" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable.HiddenDomains
done

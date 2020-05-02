#!/bin/sh
# Jake Yeung
# 2-bam_to_bigwig_ZF.sh
#  
# 2020-04-14


jmem='16G'
jtime='6:00:00'

mapq="40"

bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.NormByGenomeContentZF.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.imputevarfilt.lessstringent.mapq_${mapq}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-14.corrfilt/correlated_bins.chromofixed.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

while [[ `qstat | grep split | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

bsizes="1000 10000"
for bsize in $bsizes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bigwigs_by_cluster.mapq_${mapq}.bsize_${bsize}.NormByRPGC"
    [[ ! -d $outdir ]] && mkdir $outdir
    for inbam in `ls -d $indir/*.bam`; do
        echo $inbam
        bname=$(basename $inbam)
        bname=${bname%.*}  # strip extension 
        outbw=${outdir}/${bname}.${bsize}.bw

        [[ -e $outbw ]] && echo "$outbw found, continuing" && continue

        BNAME=$outdir/$bname.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N "bam2bw_${bname}"
    done
done

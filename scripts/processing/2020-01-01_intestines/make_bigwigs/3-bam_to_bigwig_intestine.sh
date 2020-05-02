#!/bin/sh
# Jake Yeung
# 2-bam_to_bigwig.sh
#  
# 2020-01-10

# WRAP UP
while [[ `qstat | grep split | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='6:00:00'

bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_redo_2019-12-13.split_bams_by_clusters"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_redo_2019-12-13.cluster_analysis/bams_split_by_cluster"
# indir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/raw_data/HVG-intestines.tagged_bams.by_cluster"
indir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/raw_data/HVG-intestines.tagged_bams.by_cluster"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/raw_data/HVG-intestines.bigwigs.by_cluster"
[[ ! -d $outdir ]] && mkdir $outdir
# bl="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_B6.2019-12-16/blacklist_corr_bins_merged.merged_sorted.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/from_rstudioiserver/quality_controls_intestines.OK/blacklist_corr_bins_merged.merged_sorted.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

bsize=1000
for inbam in `ls -d $indir/*.bam`; do
    echo $inbam
    bname=$(basename $inbam)
    bname=${bname%.*}  # strip extension 
    outbw=${outdir}/${bname}.bw

    [[ -e $outbw ]] && echo "$outbw found, continuing" && continue

    BNAME=$outdir/$bname.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N "bam2bw_${bname}"
done

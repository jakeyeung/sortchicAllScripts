#!/bin/sh
# Jake Yeung
# 4-make_bigwigs_from_bams.sh
#  
# 2020-10-23

jmem='16G'
jtime='2:00:00'

jbatch="batch2"
binsize=100
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file"
indir="${inmain}/split_by_cluster.MAPQ_40.${jbatch}/merged_bam"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir="${inmain}/bigwigs_BM_${jbatch}_binsize_${binsize}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/bigwigs_BM_round1_binsize_${binsize}"
[[ ! -d $outdir ]] && mkdir $outdir
bs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic-functions/scripts/processing_scripts/bam_to_bigwig.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outbw=$outdir/${bname}.bw
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $binsize"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done


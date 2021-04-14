#!/bin/sh
# Jake Yeung
# 4-make_bigwigs_from_bams.sh
#  
# 2020-10-23

jmem='16G'
jtime='2:00:00'

binsize=100
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.second_round/merged_bams"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/split_by_cluster.MAPQ_40/merged_bams"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/split_by_cluster.MAPQ_40/merged_bigwigs_binsize_${binsize}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/bigwigs_BM_round2_binsize_${binsize}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/bigwigs_BM_round1_binsize_${binsize}"
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


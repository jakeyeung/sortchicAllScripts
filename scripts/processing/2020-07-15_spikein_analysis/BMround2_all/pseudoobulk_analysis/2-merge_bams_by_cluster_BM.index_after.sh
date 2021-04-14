#!/bin/sh
# Jake Yeung
# 2-merge_bams_by_cluster_BM.sh
#  
# 2020-10-31

jmem='4G'
jtime='2:00:00'

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/split_by_cluster.MAPQ_40/merged_bams"

for b in `ls -d $indir/*.bam`; do
    bname=$(basename $b)
    bname=${bname%.*}
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools index $b"
    BNAME=${indir}/${bname}.index_sbatch_log
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${clust} --wrap "$cmd"
done


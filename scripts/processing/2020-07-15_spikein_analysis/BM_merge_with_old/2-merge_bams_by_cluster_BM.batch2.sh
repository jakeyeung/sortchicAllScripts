#!/bin/sh
# Jake Yeung
# 2-merge_bams_by_cluster_BM.sh
#  
# 2020-10-31

jmem='32G'
jtime='12:00:00'

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/split_by_cluster.MAPQ_40"
inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/split_by_cluster.MAPQ_40.batch2"
outbase="${inbase}/merged_bam"
# outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/split_by_cluster.MAPQ_40/merged_bams"
[[ ! -d $outbase ]] && mkdir $outbase

for jmark in $jmarks; do
    indir=${inbase}/${jmark}
    cd $indir
    clusts=`ls -d *.bam | cut -d"." -f4 | sort | uniq`
    echo $jmark
    echo $clusts

    for clust in $clusts; do
        BNAME=${outbase}/${jmark}_${clust}_log
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        inbams=`ls -d $indir/*${clust}.sorted.bam | tr '\n' ' '`
        outbam=${outbase}/BM_round2_all.${jmark}.${clust}.sorted.bam
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outbam $inbams; samtools index $outbam"
        echo $cmd
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${clust} --wrap "$cmd"
    done
done

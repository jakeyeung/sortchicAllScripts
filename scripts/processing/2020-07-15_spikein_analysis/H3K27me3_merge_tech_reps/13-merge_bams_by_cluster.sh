#!/bin/sh
# Jake Yeung
# 13-merge_bams_by_cluster.sh
# Merge before calling peaks 
# 2020-12-22

jmem='16G'
jtime='2:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster"
outdir="${indir}/bams_merged_by_cluster"
[[ ! -d ${outdir} ]] && mkdir ${outdir}

# # Merge Basophils and Basophils2 merge them together
# clusterscheck=`ls -d $indir/*.bam | cut -d"." -f4 | sort | uniq`
# echo $clusterscheck

# clusters="Basophils Basophils2 Bcells DCs Eryths Granulocytes HSPCs NKs pDCs"
clusters="Basophils Bcells DCs Eryths Granulocytes HSPCs NKs pDCs"

cd $indir

jmark="H3K27me3"

for cluster in $clusters; do

    BNAME=${outdir}/merge_${cluster}_sbatch
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo $cluster
    inbams=`ls -d PZ-BM-*.sorted.tagged.${cluster}*.sorted.bam | tr '\n' '\r\n'`
    ls -d PZ-BM-*.sorted.tagged.${cluster}*.sorted.bam
    # echo $check

    outbam="$outdir/PZ-BM-${jmark}-${cluster}-merged.bam"
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outbam $inbams; samtools index $outbam"
    echo $cmd
    # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=merge_$cluster --wrap "$cmd"
done


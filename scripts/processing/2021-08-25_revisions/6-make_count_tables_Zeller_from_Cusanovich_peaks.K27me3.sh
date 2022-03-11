#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=30

# # bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file"
# bamdir="${inmain}/merged_bams.first_and_second_rounds"
# outdir="${inmain}/count_tables_from_peaks"
# # outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_peaks"
# [[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams"
outdir="${inmain}/count_tables_from_Cusanovich_atac_peaks_H3K27me3"
[[ ! -d $outdir ]] && mkdir $outdir

# bedfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames_filt.mm10.nochromo.bed"
bedfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames_filt.mm10.nochromo.fourthcol.bed"
[[ ! -e $bedfile ]] && echo "HDout $bedfile not found, exiting" && exit 1
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nochr.bed"
[[ ! -e ${bl} ]] && echo "BL ${bl} not found, exiting" && exit 1

for inbam in `ls -d ${inmain}/*.bam`; do
    bambase=$(basename $inbam)
    bambase=${bambase%.*}
    bname="merged.${mark}.cutoff_analysis.merged.nochr"
    outftab=${outdir}/$bambase.${bname}.count_table.CusanovichPeaks.txt

    BNAME=${outdir}/${bname}_qsub_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $bedfile -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MakeTable_${bname}
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $bedfile -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${bname} --wrap "$cmd"
done


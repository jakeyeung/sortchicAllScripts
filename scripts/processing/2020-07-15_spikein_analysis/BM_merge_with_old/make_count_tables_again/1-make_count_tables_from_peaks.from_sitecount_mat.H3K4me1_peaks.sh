#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=40

markref="H3K4me1"

# # bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file"
# bamdir="${inmain}/merged_bams.first_and_second_rounds"
# outdir="${inmain}/count_tables_from_peaks"
# # outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_peaks"
# [[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file"
bamdir="${inmain}/merged_bams.first_and_second_rounds"
outdir="${inmain}/count_tables_from_peaks.from_sitecount_mat.${markref}_peaks"
[[ ! -d $outdir ]] && mkdir $outdir

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nochr.bed"
[[ ! -e ${bl} ]] && echo "BL ${bl} not found, exiting" && exit 1

# hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
# hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/beds.BMAllMerged2.from_peaks.fromNmat"
# hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/beds.BMAllMerged2.from_peaks.from_sitecount_mat"
# [[ ! -e $hdbase ]] && echo "HDbase $hdbase not found, exiting" && exit 1

hddir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/beds.BMAllMerged2.from_peaks.from_sitecount_mat"
[[ ! -d $hddir ]] && echo "$hddir not found, exiting" && exit 1

marks="H3K4me1 H3K27me3 H3K4me3 H3K9me3"

# hd_merged.H3K4me1.FromR.maxcount_40_60_R/
# merged.H3K4me1.cutoff_analysis.merged.nochr.bed
for mark in ${marks}; do
    for inbam in `ls -d ${bamdir}/BM_round1_round2_merged_${mark}_*.bam`; do
        bambase=$(basename $inbam)
        bambase=${bambase%.*}
        # bname="merged.${mark}.minlength_2500.cutoff_analysis.merged.nochr"
        # bname="merged.${mark}.cutoff_analysis.merged.nochr"
        bname="beds_from_sitecount_matrix.${markref}"
        # hddir="${hdbase}/hd_merged.${mark}.FromR.maxcount_40_60_R"
        bedfile="${hddir}/${bname}.bed"
        [[ ! -e $bedfile ]] && echo "HDout $bedfile not found, exiting" && exit 1
        outftab=${outdir}/$bambase.${bname}.count_table.txt

        BNAME=${outdir}/${bname}_qsub_log
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $bedfile -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MakeTable_${bname}
    done
done


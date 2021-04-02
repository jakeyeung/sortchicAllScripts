#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=40

# # bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file"
# bamdir="${inmain}/merged_bams.first_and_second_rounds"
# outdir="${inmain}/count_tables_from_peaks"
# # outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_peaks"
# [[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file"
bamdir="${inmain}/merged_bams.first_and_second_rounds"
outdir="${inmain}/count_tables_from_peaks_by_celltype"
[[ ! -d $outdir ]] && mkdir $outdir

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nochr.bed"
[[ ! -e ${bl} ]] && echo "BL ${bl} not found, exiting" && exit 1

# hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"

# hd_merged.H3K4me1.FromR.maxcount_40_60_R/
# merged.H3K4me1.cutoff_analysis.merged.nochr.bed

hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_R"

for mark in ${marks}; do
    for hddir in `ls -d $hdbase/*$mark*.cutoff`; do
        hddirbase=$(basename $hddir)  # maxcount_Rx2500.BM_round1_round2_merged_H3K27me3_DCs.2500.cutoff/
        ctype=$(echo $hddirbase | cut -d"." -f2 | cut -d"_" -f6)
        mlength=$(echo $hddirbase | cut -d"." -f3)
        echo "$mark $ctype $mlength"
        echo $bedfile
        bedfile="${hddir}/BM_round1_round2_merged_${mark}_${ctype}.${mlength}.cutoff_analysis.bed"
        [[ ! -e $bedfile ]] && echo "$bedfile not found, exiting" && exit 1
        # H3K27me3_DCs.2500.cutoff_analysis.bed"

        for inbam in `ls -d ${bamdir}/BM_round1_round2_merged_${mark}_*.bam`; do
            bambase=$(basename $inbam)
            bambase=${bambase%.*}
            outftab=${outdir}/$bambase.${mark}.${ctype}.${mlength}.count_table.txt

            BNAME=${outdir}/make_count_tables.${mark}.${ctype}.${mlength}.qsub_log
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

            # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $bedfile -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MakeTable_${bname}
            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $bedfile -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
            # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MakeTab_${ctype}.${mark}.${mlength} --wrap "$cmd"
            # echo $cmd
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MakeTab_${ctype}.${mark}.${mlength} --wrap "$cmd"
        done
    done
    # # [[ ! -e $hdbase ]] && echo "HDbase $hdbase not found, exiting" && exit 1
done


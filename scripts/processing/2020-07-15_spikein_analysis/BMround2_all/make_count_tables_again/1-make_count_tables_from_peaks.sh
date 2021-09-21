#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=40

bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_peaks"
[[ ! -d $outdir ]] && mkdir $outdir

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nochr.bed"
[[ ! -e ${bl} ]] && echo "BL ${bl} not found, exiting" && exit 1

hdbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
[[ ! -e $hdbase ]] && echo "HDbase $hdbase not found, exiting" && exit 1

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"

for mark in ${marks}; do
    for inbam in `ls -d ${bamdir}/BM_round1_round2_merged_${mark}_*.bam`; do
        bambase=$(basename $inbam)
        bambase=${bambase%.*}
        bname="merged.${mark}.minlength_2500.cutoff_analysis.merged.nochr"
        hddir="${hdbase}/hd_merged.${mark}.minlength_2500.FromR.maxcount_40_60_80"
        bedfile="${hddir}/${bname}.bed"
        [[ ! -e $bedfile ]] && echo "HDout $bedfile not found, exiting" && exit 1
        outftab=${outdir}/$bambase.${bname}.count_table.txt

        BNAME=${outdir}/${bname}_qsub_log
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $bedfile -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MakeTable_${bname}
    done
done


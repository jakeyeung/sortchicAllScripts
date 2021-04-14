#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=40

# marks="H3K9me3-H3K4me1"
marks="H3K4me1 H3K9me3"

# bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_TSS"
# [[ ! -d $outdir ]] && mkdir $outdir

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed"
# bamdir="${inmain}/tagged_bams"

bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
[[ ! -d $bamdir ]] && echo "$bamdir not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nochr.bed"
[[ ! -e ${bl} ]] && echo "BL ${bl} not found, exiting" && exit 1

# dist="1000"
# dist="10000"
# bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.nochromo.bed"
# bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.nochromo.bed"
# bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.${dist}.again.nochromo.sorted.bed"
# bedfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/beds.BMAllMerged2.from_peaks.from_sitecount_mat/beds_from_sitecount_matrix.H3K4me1.bed"
bedfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/regions_H3K4me1_H3K9me3_dynamic_regions/H3K4me1_H3K9me3_celltype_specific_genes_and_bins.2021-01-29.txt"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/count_tables_from_H3K4me1_H3K9me3_dynamic_genes_and_bins"
[[ ! -d $outdir ]] && mkdir -p $outdir

for mark in ${marks}; do
    for inbam in `ls -d ${bamdir}/BM_round1_round2_merged_*${mark}*.bam`; do
        bname=$(basename $inbam)
        bname=${bname}
        outftab=${outdir}/${bname}.count_table_k4_k9_dynamic_regions.txt
        [[ -e $outftab ]] && echo "$outftab found, continuing" && continue
        BNAME=${outdir}/${bname}_qsub_log
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed ${bedfile} -blacklist $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MakeTable_${bname}
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed ${bedfile} -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=TSS_${bname} --wrap "$cmd"
    done
done


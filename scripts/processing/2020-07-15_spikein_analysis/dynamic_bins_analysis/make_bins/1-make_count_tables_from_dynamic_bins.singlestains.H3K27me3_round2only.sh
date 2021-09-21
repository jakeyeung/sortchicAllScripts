#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=40

# marks="H3K9me3-H3K4me1"
# marks="H3K4me1 H3K4me3 H3K9me3"  # H3K27me3 done separately
marks="H3K27me3"  # H3K27me3 done separately

bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/rep2_rep3reseq_bams_together"
[[ ! -d $bamdir ]] && echo "$bamdir not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nochr.bed"
[[ ! -e ${bl} ]] && echo "BL ${bl} not found, exiting" && exit 1

beddir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.dynamic_bins_TSS_TES_regions"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/count_tables_all_four_marks_dynamic_bins"
[[ ! -d $outdir ]] && mkdir -p $outdir

bedfnames="dynamic_bins.50kb.H3K27me3.txt celltype_specific_genes.TSS_TES.txt celltype_specific_genes.TSS_10kb.txt"

for mark in ${marks}; do
    for bedfname in $bedfnames; do
        bedfile=${beddir}/${bedfname}
        echo $bedfile
        [[ ! -e $bedfile ]] && echo "$bedfile not found, exiting" && exit 1
        bedname=$(basename $bedfile)
        bedname=${bedname%.*}

        for inbam in `ls -d ${bamdir}/PZ*${mark}*.bam`; do
            bname=$(basename $inbam)
            bname=${bname}
            outftab=${outdir}/${bname}.${bedname}.txt
            [[ -e $outftab ]] && echo "$outftab found, continuing" && continue
            BNAME=${outdir}/${bname}_qsub_log
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed ${bedfile} -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=TSS_${bname} --wrap "$cmd"
        done
    done
done


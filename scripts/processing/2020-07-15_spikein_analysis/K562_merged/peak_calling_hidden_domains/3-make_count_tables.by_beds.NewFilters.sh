#!/bin/sh
# Jake Yeung
# 6-make_count_tables.sh
# Make count tables
# 2020-01-02


# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_CellCycleSorted_split_by_cellcycle.NoChrPrefix"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="${inmain}/countTablesAndRZr1only_ByBed.NewFilters"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='16G'
jtime='6:00:00'

mapq=40
# binsize=100000
# stepsize=$binsize
# stepsize=50000

# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
# inbed="/hpc/hub_oudenaarden/avo/scChiC/published_ChIP_data/H3K27me3_ENCFF000VDN_unique_sorted.deduplicated.bed"
# inbed="/hpc/hub_oudenaarden/jyeung/data/public_data/H3K27me3_ENCFF000VDN_unique_sorted.deduplicated.labeled_peaks.bed"
# inbedmain="/hpc/hub_oudenaarden/jyeung/data/public_data/ENCODE/bedfiles"

# K562_AllMerged_H3K27me3.merged.sorted.tagged.G1filt.sorted.1000.cutoff_analysis.bed"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

cellcycles="X0_G1 X1_S X2_G2_M"

for jmark in $jmarks; do
    inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt_NoChrPrefix/hiddendomains_outputs/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.1000.cutoff"
    for cc in $cellcycles; do
        # inbed="${inbedmain}/ENCODEpeaks.${jmark}.nochr.cleanchr.qval_${qvalcutoff}.bed"
        inbed="${inbedmain}/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.1000.cutoff_analysis.bed"
        # bname="K562_CellCycleSorted_${jmark}.merged.sorted.tagged.bam"
        bname="K562_CellCycleSorted_${jmark}.merged.sorted.tagged.${cc}.sorted.bam"
        inbam=${inmain}/${bname}
        [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1

        outf=$outdir/${bname}.countTable.bedfile.csv
        [[ -e $outf ]] && echo "$outf found, continuing" && continue

        BNAME=$outdir/${bname}.counttables.ByBed.again2.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $inbed"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts_${bname} --wrap "$cmd"
    done
done


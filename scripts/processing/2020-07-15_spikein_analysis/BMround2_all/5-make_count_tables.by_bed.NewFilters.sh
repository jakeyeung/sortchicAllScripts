#!/bin/sh
# Jake Yeung
# 6-make_count_tables.sh
# Make count tables
# 2020-01-02

# WRAP UP
while [[ `squeue -u jyeung | grep merge_ | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done


# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/merged_bams"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="${inmain}/countTablesAndRZr1only_ByChromo.NewFilters.bybed.blfix"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='8G'
jtime='2:00:00'

mapq=40
# stepsize=50000
# binsize=50000

# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.nochr.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1
# inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl/good_bins_dbl.bed"
inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all.dbl.blfix/good_bins_dbl.blfix.nochromo.bed"

for inbam in `ls -d $inmain/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outf=$outdir/${bname}.countTable.dbl_bins.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${bname}.counttables.ByChromo.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $inbed"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

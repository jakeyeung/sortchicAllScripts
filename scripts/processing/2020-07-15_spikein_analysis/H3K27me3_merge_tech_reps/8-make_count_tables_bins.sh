#!/bin/sh
# Jake Yeung
# 8-make_count_tables_bins.sh
#  
# 2020-11-28

jmem='16G'
jtime='12:00:00'

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams"

outdir="${inmain}/counts_tables"
[[ ! -d $outdir ]] && mkdir $outdir

# jmark="H3K27me3"

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"

binsize=50000
stepsize=$binsize
# echo $jmark

for inbam in `ls -d $inmain/*.bam`; do
    inbambase=$(basename $inbam)
    inbambase=${inbambase%.*}
    echo $inbambase
    bname=${inbambase}

    BNAME=$outdir/${bname}.sbatchout
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "dbase $DBASE not found, exiting" && exit 1

    outf1=$outdir/${bname}.countTable.binsize_${binsize}.csv
    [[ -e $outf1 ]] && echo "outf1 $outf1 found, continuing" && continue

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam -sliding $stepsize --filterXA -minMQ $mapq -o $outf1 -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

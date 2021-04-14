#!/bin/sh
# Jake Yeung
# 6-make_count_tables.sh
# Make count tables
# 2020-01-02


# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="${inmain}/TSS_count_tables.NewFilters"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='32G'
jtime='6:00:00'

mapq=40
# bsize="50000"
bsize="10000"
# bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.nochromo.bed"
bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.${bsize}.again.nochromo.sorted.bed"

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

for inbam in `ls -d $inmain/*.tagged.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outf=$outdir/${bname}.countTable.TSS.bsize_${bsize}.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${bname}.counttables.ByChromo.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam  --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $bedfile"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

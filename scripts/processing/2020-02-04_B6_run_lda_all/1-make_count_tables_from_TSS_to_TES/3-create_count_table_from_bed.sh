#!/bin/sh
# Jake Yeung
# 2-create_count_table_from_bed.sh
#  
# 2020-08-10

# inf="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/5-make_count_tables.NewFilters.sh"

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir="${indir}/countTablesAndRZr1only.NewFilters"
[[ ! -d $outdir ]] && mkdir $outdir

# binsize=50000
# stepsize=50000
mapq=40

jmem='64G'
jtime='12:00:00'

# inbed="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.merged.bed"
inbed="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.merged.rearranged.bed"

for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    BNAME=${outdir}/${bname}.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf="${outdir}/${bname}.count_table.TSS_TES.csv"
    [[ ! -e $outf ]] && echo "$outf not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bedfile $inbed -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=count_${bname} --wrap "$cmd"
done


#!/bin/sh
# Jake Yeung
# 2-get_reads_at_CTCF_motifs.sh
#  
# 2020-08-17

jmem='16G'
jtime='3:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/tagged"
outdir="${indir}/count_tables_CTCF"
[[ ! -d $outdir ]] && mkdir $outdir

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/beds/CTCF_motifs/mm10_bonemarrow-all_CTCF-chip_optimal-IDRpeaks_ENCFF806PDR.chromofilt.bed"
# inbed="/hpc/hub_oudenaarden/jyeung/data/databases/beds/CTCF_motifs/mm10_bonemarrow-all_CTCF-chip_optimal-IDRpeaks_ENCFF806PDR.chromofilt.bed"
# outbed="/hpc/hub_oudenaarden/jyeung/data/databases/beds/CTCF_motifs/mm10_bonemarrow-all_CTCF-chip_optimal-IDRpeaks_ENCFF806PDR.chromofilt.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"

mapq=40

for inbam in `ls -d $indir/*.bam`; do

    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchout
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.count_tables_CTCF.csv

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $inbed"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

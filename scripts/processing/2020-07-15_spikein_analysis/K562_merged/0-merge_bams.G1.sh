#!/bin/sh
# Jake Yeung
# 0-merge_bams.sh
#  
# 2020-09-14

jmem='16G'
jtime='3:00:00'

# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams"
indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/merged_bams"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    # b1="K562-EtOH-${jmark}-G1-G2.sorted.tagged.bam"
    b1="K562-EtOH-${jmark}.merged.sorted.tagged.bam"
    b2="K562-EtOH-${jmark}.CellCycleSorted.merged.sorted.tagged.bam"
    b3="K562-EtOH-${jmark}.G1sorted.merged.sorted.tagged.bam"

    inf1=${indir1}/${b1}
    inf2=${indir2}/${b2}
    inf3=${indir2}/${b3}

    [[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1
    [[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1
    [[ ! -e $inf3 ]] && echo "$inf3 not found, exiting" && exit 1

    outf="${outdir}/K562_AllMerged_${jmark}.merged.sorted.tagged.bam"

    BNAME="${outdir}/K562_AllMerged_${jmark}.log"
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $inf1 $inf2"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=merge_CellCycle_${jmark} --wrap "$cmd"
done

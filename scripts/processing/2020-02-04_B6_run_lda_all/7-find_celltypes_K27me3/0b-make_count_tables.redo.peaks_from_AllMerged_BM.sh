#!/bin/sh
# Jake Yeung
# 0b-make_count_tables.sh
#  
# 2020-04-03

jmem='16G'
jtime='3:00:00'
ncores=1

mapq=40

jmarkref="H3K4me1"
# jmarkref="H3K27me3"
inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/bed_from_count_mat/${jmarkref}_peaks_from_hiddendomains.nochromo.bed"
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

blacklist="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.nochr.bed"

indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_mflorescu/B6_unfixed"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix"
outdirtab="$outdir/countTables.unfixed"
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $outdirtab ]] && mkdir $outdirtab

jmarks="K27m3 K4m1_K27m3 K4m1"

for jmark in $jmarks; do
    inbam="$indir/all_BM_${jmark}_200119.bam"
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    bname=$(basename $inbam)
    bname=${bname%.*}
    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outftab=$outdirtab/$bname.mq_${mapq}.on_${jmarkref}_peaks.csv
    jobname=${bname}

    cmdcounts=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name -binTag DS --dedup -blacklist $blacklist --r1only -bedfile $inbed; gzip $outftab"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=$ncores --nodes=1 --ntasks-per-node=$ncores --ntasks-per-socket=$ncores --job-name=${jmark}.$jobname --wrap "$cmdcounts"
done


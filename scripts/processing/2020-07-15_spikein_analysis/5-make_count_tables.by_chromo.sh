#!/bin/sh
# Jake Yeung
# 6-make_count_tables.sh
# Make count tables
# 2020-01-02


# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562.tagged_bams_mergedbymarks"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed"
outdir="${inmain}/countTablesAndRZr1only_ByChromo"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='32G'
jtime='6:00:00'

mapq=40
# stepsize=50000
# binsize=50000

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.nochr.bed"

for inbam in `ls -d $inmain/*.tagged.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outf=$outdir/${bname}.countTable.ByChromo.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${bname}.counttables.ByChromo.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup --r1only -blacklist $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup --r1only -blacklist $bl"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

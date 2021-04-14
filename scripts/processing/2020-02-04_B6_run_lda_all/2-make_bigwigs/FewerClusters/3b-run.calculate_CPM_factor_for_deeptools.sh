#!/bin/sh
# Jake Yeung
# 3b-run.calculate_default_norm_factor_for_deeptools.sh
#  
# 2020-06-19

jmem='8G'
jtime='0:30:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/2-make_bigwigs/FewerClusters/calculate_CPM_factor_for_deeptools.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits.Downsamp"

# outdir="${indir}/default_norm_factors"
outdir="${indir}/CPM_factors"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits.Downsamp/default_norm_factors"
[[ ! -d $outdir ]] && mkdir $outdir

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1


for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatch
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.DefaultNormFactor.txt
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -bam $inbam -outfile $outf -blackListFileName $bl"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -bam $inbam -outfile $outf -blackListFileName $bl
    # echo "Done testing"
    # exit 0
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
done



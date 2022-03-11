#!/bin/sh
# Jake Yeung
# 3-make_bigwigs_from_bams.sh
#  
# 2022-01-10

# WRAP UP
while [[ `squeue -u jyeung | grep K562_rev | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/K562/merged_bams/bams_split_by_good_cells_NoChrPrefix"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/K562/bigwigs"

bs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic-functions/scripts/processing_scripts/bam_to_bigwig.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

binsize=1000
for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outbw=$outdir/${bname}.bw
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $binsize"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done


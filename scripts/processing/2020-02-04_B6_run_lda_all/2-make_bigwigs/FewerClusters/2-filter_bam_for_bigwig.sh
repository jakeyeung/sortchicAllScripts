#!/bin/sh
# Jake Yeung
# 2-filter_bam_for_bigwig.sh
# Remove duplicates, MAPQ40, R1only, blacklist, not atlernative hits
# 2020-06-18

jmem='16G'
jtime='2:00:00'


mapq=40
ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-06-03/singlecellmultiomics/bamProcessing/bamFilter.py"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_${mapq}.2020-06-17"
outdir="$indir/DedupR1onlyNoAltHits"

[[ ! -d $outdir ]] && mkdir $outdir

exprs='r.is_read1 and not r.is_duplicate and not read_has_alternative_hits_to_non_alts(r)'

for inbam in `ls -d $indir/*bam`; do
    bname=$(basename $inbam)
    bnamestripped=${bname%.*}

    BNAME=${outdir}/${bnamestripped}.sbatch
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outbam=${outdir}/${bnamestripped}.cleaned.bam
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam \"$exprs\" -o $outbam; samtools index $outbam"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bnamestripped} --wrap "$cmd"

    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
    # echo "python $ps $inbam \"$exprs\" -o $outbam"
    # exit 0

done



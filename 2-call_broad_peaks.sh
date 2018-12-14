#!/bin/sh
# Jake Yeung
# 2-call_broad_peaks.sh
# Call broad peaks  
# 2018-12-10

# Call broad peaks

jmem='10G'
jtime='3:00:00'

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods"
badbam="/hpc/hub_oudenaarden/jyeung/data/histone-mods/PZ-BM-m1-H3K9me3-2_AH3VGVBGX9_S1/PZ-BM-m1-H3K9me3-2_AH3VGVBGX9_S1.filtered.sorted.bam"  # an empty bam file skip it
outdir=$inmain

for indir in $(ls -d $inmain/PZ-BM*-H*H*_S*); do
    # indir="/hpc/hub_oudenaarden/jyeung/data/histone-mods/PZ-BM-m1-H3K4me3-1_AH3VGVBGX9_S10"
    bname=$(basename $indir)
    inf=$indir/$bname.filtered.sorted.bam
    [[ "$inf" == "$badbam" ]] && echo "$inf is bad, skipping" && continue
    [[ ! -e $inf ]] && echo "$inf not found, stopping" && exit 1
    outdir="$indir/macs2_out"
    [[ ! -d $outdir ]] && mkdir $outdir
    bname=$(basename $inf)
    bname=${bname%%.*}
    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo "conda activate py2; macs2 callpeak -t $inf --broad -g mm --broad-cutoff 0.1 --min-length 1000 -n $bname --outdir $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
     # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2; macs2 callpeak -t $inf --broad -g mm --pvalue 0.01 --broad-cutoff 0.01 --min-length 500 -n $bname --outdir $outdir" 
     # echo ${BNAME}.out
     # echo ${BNAME}.err
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2; macs2 callpeak -t $inf --broad -g mm --pvalue 0.01 --broad-cutoff 0.01 --min-length 500 -n $bname --outdir $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
done

#!/bin/sh
# Jake Yeung
# 7-merge_bams_by_marks.K4me1.sh
#  
# 2020-04-07

jmem='32G'
jtime='12:00:00'

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
cd $indir

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks"
[[ ! -d $outdir ]] && mkdir $outdir

for mark in $marks; do
    # mark="H3K4me1"
    bin=`ls $indir/*.bam | grep -i ${mark} | tr "\n" " "`
    bout="$outdir/PZ-ChIC-ZF_${mark}_2020-04-07.bam"

    BNAME="$outdir/merge_${mark}_qsub"
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo $bin
    echo $bout
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $bout $bin --output-fmt BAM; samtools index $bout" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N merge_${mark}
done


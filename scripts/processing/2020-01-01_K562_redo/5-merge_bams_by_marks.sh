#!/bin/sh
# Jake Yeung
# 5-merge_bams_by_marks.sh
# Merge wildl types together. 
# 2019-12-25

jmem='96G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562/PZ-K562.tagged_bams"
cd $indir
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562.tagged_bams_mergedbymarks"

# merge H3K4me1, remove rep3 
# rep3 was renamed 
# mv PZ-K562-G1-H3K4me1-rep3.tagged.bam PZ-K562-G1-H3K4me1-BADREP-rep3.tagged.bam
infs=$indir/PZ-K562-G1-H3K4me1-rep*.tagged.bam
outf=$outdir/PZ-K562-G1-H3K4me1-Rep1Rep2merged.tagged.bam
echo $infs
echo $outfs
BNAME=$outdir/H3K4me1.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

# # merge H3K4me3, all reps
infs=$indir/PZ-K562-G1-H3K4me3-rep*.tagged.bam
outf=$outdir/PZ-K562-G1-H3K4me3-merged.tagged.bam
echo $infs
echo $outfs
BNAME=$outdir/H3K4me3.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1


# merge H3K27me3, all reps
infs=$indir/PZ-K562-G1-H3K27me3-rep*.tagged.bam
outf=$outdir/PZ-K562-G1-H3K27me3-merged.tagged.bam
echo $infs
echo $outfs
BNAME=$outdir/H3K27me3.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

# merge H3K9me3, all reps 
infs=$indir/PZ-K562-G1-H3K9me3-rep*.tagged.bam
outf=$outdir/PZ-K562-G1-H3K9me3-merged.tagged.bam
echo $infs
echo $outfs
BNAME=$outdir/H3K9me3.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outf $infs; samtools index $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

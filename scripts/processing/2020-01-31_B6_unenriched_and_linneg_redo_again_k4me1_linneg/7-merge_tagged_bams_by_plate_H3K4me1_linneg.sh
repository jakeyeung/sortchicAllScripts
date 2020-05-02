#!/bin/sh
# Jake Yeung
# 7-merge_tagged_bams_by_plate_H3K4me1.sh
# 
# 2019-12-14



indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_linneg_K4me1.tagged_bams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_linneg_K4me1.tagged_bams_mergedbymarks"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='64G'
jtime='6:00:00'
BNAME=$outdir/mergebams_qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

# Merge H3K4me1: Unenriched and StemCells (lineage neg not yet done)
bin1=/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/H3K4me1-BM_SC-merged.tagged.bam
bin2=$indir/*H3K4me1*.retagged.bam
bin2="$indir/PZ-ChIC-Bl6-BM-lin-H3K4me1-1.retagged.bam $indir/PZ-ChIC-Bl6-BM-lin-H3K4me1-2.retagged.bam $indir/PZ-ChIC-Bl6-BM-lin-H3K4me1-3.retagged.bam"
for b in $bin2; do
    echo $b; 
    [[ ! -e $b ]] && echo "$b not found, exiting" && exit 1
done
# bin="$bin1 $bin2"
bout="$outdir/H3K4me1-BM_Linneg_SC-merged.2020-01-31.tagged.bam"

echo "samtools merge $bout $bin1 $bin2 --output-fmt BAM; samtools index $bout" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N mergebams


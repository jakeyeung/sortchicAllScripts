#!/bin/sh
# Jake Yeung
# 0-normalize_bam_to_input.sh
#  
# 2020-11-22

jmem='16G'
jtime='2:00:00'

indir="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/chipseq_from_ENCODE_first_submission"

binput="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/input_unique_sorted.deduplicated.bam"
b1="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/H3K4me1_ENCFF063EDR.deduplicated.bam"
b2="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/H3K4me3_ENCFF000VDX_unique_sorted.deduplicated.bam"
b3="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/H3K27me3_ENCFF000VDN_unique_sorted.deduplicated.bam"
b4="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/K9me3_unique_sorted.deduplicated.bam"
bsize=1000

BNAME=${outdir}/H3K4me1_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b1 -b2 $binput -o $outdir/H3K4me1_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K4me3_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b2 -b2 $binput -o $outdir/H3K4me3_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K27me3_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b3 -b2 $binput -o $outdir/H3K27me3_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K9me3_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b4 -b2 $binput -o $outdir/H3K9me3_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

#!/bin/sh
# Jake Yeung
# 0-normalize_bam_to_input.sh
#  
# 2020-11-22

jmem='16G'
jtime='2:00:00'

# indir="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/chipseq_from_ENCODE_newbams_norm_by_input"
[[ ! -d $outdir ]] && mkdir $outdir

# binput="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/input_unique_sorted.deduplicated.bam"
binput="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/input/ENCFF994FIB_unique_sorted.bam"
# b1="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/H3K4me1_ENCFF063EDR.deduplicated.bam"
b1a="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K4me1/ENCFF000BYG_unique_sorted.bam"
b1b="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K4me1/ENCFF000BXX_unique_sorted.bam"


# b2="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/H3K4me3_ENCFF000VDX_unique_sorted.deduplicated.bam"
b2a="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K4me3/ENCFF894KBP_unique_sorted.bam"
b2b="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K4me3/ENCFF010SAE_unique_sorted.bam"

# b3="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/H3K27me3_ENCFF000VDN_unique_sorted.deduplicated.bam"
b3a="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K27me3/ENCFF000BXP_unique_sorted.bam"
b3b="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K27me3/ENCFF000BXN_unique_sorted.bam"

# b4="/hpc/hub_oudenaarden/Peter/data/K562/published/published_for_paper/K9me3_unique_sorted.deduplicated.bam"
b4a="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K9me3/ENCFF001QWX_unique_sorted.bam"
b4b="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams/H3K9me3/ENCFF001QWW_unique_sorted.bam"

[[ ! -e $b1a ]] && echo "$b1a not found, exiting" && exit 1
[[ ! -e $b1b ]] && echo "$b1b not found, exiting" && exit 1

[[ ! -e $b2a ]] && echo "$b2a not found, exiting" && exit 1
[[ ! -e $b2b ]] && echo "$b2b not found, exiting" && exit 1

[[ ! -e $b3a ]] && echo "$b3a not found, exiting" && exit 1
[[ ! -e $b3b ]] && echo "$b3b not found, exiting" && exit 1

[[ ! -e $b4a ]] && echo "$b4a not found, exiting" && exit 1
[[ ! -e $b4b ]] && echo "$b4b not found, exiting" && exit 1

bsize=1000

BNAME=${outdir}/H3K4me1a_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b1a -b2 $binput -o $outdir/H3K4me1_a_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K4me1b_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b1b -b2 $binput -o $outdir/H3K4me1_b_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K4me3a_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b2a -b2 $binput -o $outdir/H3K4me3_a_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K4me3b_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b2b -b2 $binput -o $outdir/H3K4me3_b_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K27me3a_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b3a -b2 $binput -o $outdir/H3K27me3_a_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K27me3b_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b3b -b2 $binput -o $outdir/H3K27me3_b_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K9me3a_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b4a -b2 $binput -o $outdir/H3K9me3_a_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

BNAME=${outdir}/H3K9me3b_normalize.${bsize}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamCompare -b1 $b4b -b2 $binput -o $outdir/H3K9me3_b_input_normalized.${bsize}.bw --binSize $bsize" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1




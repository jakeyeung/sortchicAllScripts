#!/bin/sh
# Jake Yeung
# 5-run.make_count_from_bam.sh
#  
# 2021-08-26

jmem='16G'
jtime='3:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-08-25_revisions/make_count_from_bam.py"
# inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_remapped_mm10/BoneMarrow_62216.bam.remapped.sorted.bam"
inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_remapped_mm10/BoneMarrow_62016.bam.remapped.sorted.bam"

infcoords="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.10000.again.nochromo.sorted.bed"
# infcoords="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.merged_with_blacklist.nochromo.sorted.bed"
infcells="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.noheader.txt"

mapq=40
bname=$(basename $inbam)
bname=${bname%.*}
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables/"
# outfile="${outdir}/${bname}.TSS_50kb.mapq_${mapq}.count_table.txt"
outfile="${outdir}/${bname}.TSS_10kb.mapq_${mapq}.count_table.txt"


BNAME=${outdir}/${bname}.sbatch_output
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -inbam $inbam -metadata_cells $infcells -metadata_coords $infcoords -outfile $outfile -mapq $mapq"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${bname} --wrap "$cmd"

# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -inbam $inbam -metadata_cells $infcells -metadata_coords $infcoords -outfile $outfile

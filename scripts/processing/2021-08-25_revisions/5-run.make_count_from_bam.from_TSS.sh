#!/bin/sh
# Jake Yeung
# 5-run.make_count_from_bam.sh
#  
# 2021-08-26

jmem='96G'
jtime='1:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-08-25_revisions/make_count_from_bam.py"

mapq=0
bnames="BoneMarrow_62216 BoneMarrow_62016"

# inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_remapped_mm10/BoneMarrow_62016.bam.remapped.sorted.bam"
# infcoords="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.10000.again.nochromo.sorted.bed"
# infcoords="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.merged_with_blacklist.nochromo.sorted.bed"

# infcoords="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames.mm10.nochromo.bed"
# infcoords="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames_filt.mm10.nochromo.bed"
# infcoords="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.10000.again.nochromo.sorted.bed"
infcoords="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.nochromo.bed"
# infcells="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.noheader.txt"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/count_tables"

for bname in $bnames; do

    inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_remapped_mm10/${bname}.bam.remapped.sorted.bam"
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    infcells="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/metadata/cell_metadata.bonemarrow_only.noheader.${bname}.txt"
    [[ ! -e $infcells ]] && echo "$infcells not found, exiting" && exit 1

    outfile="${outdir}/${bname}.50kb_TSS.mapq_${mapq}.count_table.txt"

    BNAME=${outdir}/${bname}.sbatch_output
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -inbam $inbam -metadata_cells $infcells -metadata_coords $infcoords -outfile $outfile -mapq $mapq --quiet"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -inbam $inbam -metadata_cells $infcells -metadata_coords $infcoords -outfile $outfile -mapq $mapq --quiet
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${bname} --wrap "$cmd"
done

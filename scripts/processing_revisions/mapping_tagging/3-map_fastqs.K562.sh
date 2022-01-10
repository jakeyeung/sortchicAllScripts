#!/bin/sh
# Jake Yeung
# 3-map_fastq.sh
#  
# 2022-01-06

jmem='32G'
jtime="48:00:00"
ncores=4

# # WRAP UP
# while [[ `squeue -u jyeung | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

bwabin="bwa"
# indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_97_129B6Masked_ERCC92.fa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92.fa"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/raw_demultiplexed"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/raw_demultiplexed/K562"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for indir in `ls -d $inmain/PZ-sortChIC-k562*`; do
    bname=$(basename $indir)
    outdir=$indir

    f1=$indir/trimmed.R1.fastq.gz
    f2=$indir/trimmed.R2.fastq.gz
    [[ ! -e $f1 ]] && echo "$f1 not found, continuing" && continue
    [[ ! -e $f2 ]] && echo "$f2 not found, continuing" && continue

    BNAME=$outdir/$bname.mapping.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=$outdir/bwaMapped.bam
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=BwaMap_${bname} --wrap "$cmd"
done

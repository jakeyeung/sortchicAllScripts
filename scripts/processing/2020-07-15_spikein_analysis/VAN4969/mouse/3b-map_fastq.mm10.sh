#!/bin/sh
# Jake Yeung
# 3-map_fastq.sh
#  
# 2019-12-19

# # WRAP UP
# while [[ `squeue -u jyeung | grep jyeung |  wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

bwabin="bwa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_NOMASK_ERCC92.WithLambdaPhage.fa"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

jmem='32G'
jtime="24:00:00"
ncores=8

for indir in `ls -d $inmain/mouse*`; do
    bname=$(basename $indir)
    outdir=$indir
    f1=$indir/trimmed.R1.scmo3.fastq.gz
    f2=$indir/trimmed.R2.scmo3.fastq.gz
    [[ ! -e $f1 ]] && echo "$f1 not found, exiting" && exit 1
    [[ ! -e $f2 ]] && echo "$f2 not found, exiting" && exit 1

    BNAME=$outdir/$bname.mapping.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=$outdir/bwaMapped.bam
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=BwaMap_${bname} --wrap "$cmd"
done

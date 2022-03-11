#!/bin/sh
# Jake Yeung
# 3-map_fastq.sh
#  
# 2019-12-19

# WRAP UP
while [[ `squeue -u jyeung | grep jyeung |  wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

bwabin="bwa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92.fa"
# indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92_WithLambdaPhage.fa"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562/raw_demultiplexed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/raw_demultiplexed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/raw_demultiplexed"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs_trimmed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/bams"
[[ ! -d $outmain ]] && mkdir $outmain

jmem='32G'
jtime="24:00:00"
ncores=8

for indir in `ls -d $inmain/SRR*`; do
    bname=$(basename $indir)
    # echo $indir
    outdir=${outmain}/${bname}
    [[ ! -d $outdir ]] && mkdir $outdir
    f1="$indir/${bname}_1.trimmed.fastq.gz"
    [[ ! -e $f1 ]] && echo "$f1 not found, exiting" && exit 1

    BNAME=${outdir}/${bname}.mapping.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.bam
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f1 | samtools view -Sb - > $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=BwaMap_${bname} --wrap "$cmd"
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf"
done

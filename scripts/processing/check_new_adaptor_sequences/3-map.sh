#!/bin/sh
# Jake Yeung
# 3-map.sh
#  
# 2019-09-04

bwabin="bwa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92.fa"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/raw_demultiplexed"

jmem='32G'
jtime='12:00:00'
ncores=4

for indir in `ls -d $inmain/PZ*`; do
    bname=$(basename $indir)
    echo $indir
    outdir=$indir
    f1=$indir/trimmed.R1.fastq.gz
    f2=$indir/trimmed.R2.fastq.gz
    [[ ! -e $f1 ]] && echo "$f1 not found, exiting" && exit 1
    [[ ! -e $f2 ]] && echo "$f2 not found, exiting" && exit 1

    BNAME=$outdir/$bname.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # outf=$outdir/bwaMapped.${ncores}.bam
    outf=$outdir/bwaMapped.bam
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # bwa mem -t 4 -M ${indxbase} $f1 $f2 | samtools view -b - > ./unsorted.bam; samtools sort -T ./temp_sort -@ 4 ./unsorted.bam > ./sorted.unfinished.bam; 
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate singlecellmultiomicsenv2; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname
    # echo "echo . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate singlecellmultiomicsenv2; $bwabin mem -t 16 $indxbase $f1 $f2 | samtools view -Sb - > $outf | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname"
done

# 
# "echo "echo $outdir; $bwabin mem -t 16 $indxbase $f1 $f2 | samtools view -Sb - > $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu"
# 
# mv ./sorted.unfinished.bam ./sorted.bam; 
# samtools index ./sorted.bam; rm ./unsorted.bam;

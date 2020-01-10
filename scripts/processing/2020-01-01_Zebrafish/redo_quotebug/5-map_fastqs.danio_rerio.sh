#!/bin/sh
# Jake Yeung
# 5-map_fastqs.sh
#  
# 2019-09-28

# # WRAP UP
# while [[ `qstat | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

bwabin="bwa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/98/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/raw_demultiplexed"

jmem='32G'
jtime='12:00:00'
ncores=4

for indir in `ls -d $inmain/PZ-ChIC-ZFWKM*`; do
    bname=$(basename $indir)
    echo $indir
    outdir=$indir
    f1=$indir/trimmed.R1.fastq.gz
    f2=$indir/trimmed.R2.fastq.gz
    [[ ! -e $f1 ]] && echo "$f1 not found, exiting" && exit 1
    [[ ! -e $f2 ]] && echo "$f2 not found, exiting" && exit 1

    BNAME=$outdir/$bname.mapping.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=$outdir/bwaMapped.bam
    # [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf"
done

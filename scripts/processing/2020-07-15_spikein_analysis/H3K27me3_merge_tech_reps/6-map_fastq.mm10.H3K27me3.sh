#!/bin/sh
# Jake Yeung
# 3-map_fastq.sh
#  
# 2019-12-19

# WRAP UP
# WRAP UP
# while [[ `qstat | wc -l` > 0 ]]; do
#         echo "sleep for 600 seconds"
#         sleep 600
# done

while [[ `squeue -u jyeung | grep jyeung |  wc -l` > 1 ]]; do
        echo "sleep for 66 seconds waiting for trim"
        sleep 66
done

bwabin="bwa"
# indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92_WithLambdaPhage.fa"
indxbase="/hpc/hub_oudenaarden/group_references/ensembl/97/mus_musculus/primary_assembly_NOMASK_ERCC92.WithLambdaPhage.fa"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/raw_demultiplexed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/raw_demultiplexed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

jmem='32G'
jtime="24:00:00"
ncores=8

for indir in `ls -d $inmain/PZ*`; do
    bname=$(basename $indir)
    # echo $indir
    outdir=$indir
    # f1=$indir/trimmed.R1.fastq.gz
    # f2=$indir/trimmed.R2.fastq.gz
    f1=$indir/trimmed.R1.fastq.gz
    f2=$indir/trimmed.R2.fastq.gz
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

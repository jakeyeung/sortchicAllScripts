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
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/trimmed_r1"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='32G'
jtime="24:00:00"
ncores=8

for f1 in `ls -d $inmain/*fastq.gz`; do
    bname=$(basename $f1)
    bname=${bname%%.*}
    # echo $indir
    # outdir=${outmain}/${bname}
    # [[ ! -d $outdir ]] && mkdir $outdir
    f1="$inmain/${bname}.demux.trimmed.fastq.gz"
    [[ ! -e $f1 ]] && echo "$f1 f1 not found, exiting" && exit 1

    BNAME=${outdir}/${bname}.mapping.sbatch_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.bam
    [[ -e $outf ]] && echo "$outf outf found, continuing" && continue

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; $bwabin mem -t $ncores $indxbase $f1 | samtools view -Sb - > $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=BwaMap_${bname} --wrap "$cmd"
done

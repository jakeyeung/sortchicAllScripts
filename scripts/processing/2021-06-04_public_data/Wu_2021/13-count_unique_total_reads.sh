#!/bin/sh
# Jake Yeung
# 14-count_unique_total_reads.sh
#  
# 2021-06-15

jmem='16G'
jtime='2:00:00'

# ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/Ku_2021/count_unique_total_reads.py"
ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/Wu_2021/count_unique_total_reads_Wu.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/bams_demux_bugfix"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1_notrim_bowtie2"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/counts_output_r1_notrim_bowtie2"
[[ ! -d $outdir ]] && mkdir $outdir
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output"

for inf in `ls -d $indir/*.bam`; do
# for inf in `ls -d $indir/SRR10387113*.bam`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    outf="${outdir}/${bname}_dupcounts.bed.gz"
    BNAME="${outdir}/${bname}_dupcounts.sbatch.log"
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf -species human -mapq 0"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf -species human -mapq 0
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf
done

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/bams_demux_bugfix/SRR10615134_1.sorted.bam"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output/out.txt.gz"





#!/bin/sh
# Jake Yeung
# 10-demux_fastqs.sh
#  
# 2021-06-14

jmem='4G'
jtime='1:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/Ku_2021/demux_fastqs_Ku.py"

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs_trimmed"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs"

bcf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/GSE139857_barcode_96.txt"
[[ ! -e $bcf ]] && echo "$bcf not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs_annot_bugfix"
[[ ! -d $outdir ]] && mkdir $outdir

for d in `ls -d ${indir1}/SRR*`; do
    dbase=$(basename $d)

    BNAME=${outdir}/${dbase}.sbatch_output
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


    echo $d
    echo ${dbase}

    inf1=${indir1}/${dbase}/${dbase}_1.trimmed.fastq.gz
    inf2=${indir2}/${dbase}/${dbase}_2.fastq.gz
    [[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1
    [[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1
    outf="${outdir}/${dbase}_1.demux.trimmed.fastq.gz"
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile1 $inf1 -infile2 $inf2 -bcfile $bcf -outfile $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${dbase} --wrap "$cmd"
done

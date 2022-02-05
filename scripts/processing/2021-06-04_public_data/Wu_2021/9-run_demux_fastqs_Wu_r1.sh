#!/bin/sh
# Jake Yeung
# 9-run_demux_fastqs_Wu.sh
#  
# 2021-06-16

jmem='16G'
jtime='1:00:00'

# ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/Wu_2021/demux_fastqs_Wu.py"
ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/Wu_2021/demux_fastqs_Wu_r2.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

bname="SRR12638101"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/${bname}"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/${bname}/demux_R2"
[[ ! -d $outdir ]] && mkdir $outdir

BNAME=${outdir}/${bname}.sbatch_log
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

f1name="${bname}_1.fastq.gz"
f2name="${bname}_2.fastq.gz"

f1=${indir}/${f1name}
f2=${indir}/${f2name}

fout="${outdir}/${bname}_1.demux.fastq.gz"

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile1 $f1 -infile2 $f2 -outfile $fout"

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"

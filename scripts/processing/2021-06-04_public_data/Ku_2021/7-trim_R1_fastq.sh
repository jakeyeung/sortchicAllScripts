#!/bin/sh
# Jake Yeung
# 7-trim_R1_fastq.sh
#  
# 2021-06-09

jmem='8G'
jtime='12:00:00'

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs_trimmed"
[[ ! -d $outmain ]] && mkdir $outmain

for d in `ls -d $inmain/SRR*`; do
    dbase=$(basename $d)
    outdir=${outmain}/${dbase}
    [[ ! -d $outdir ]] && mkdir $outdir
    
    f1="${d}/${dbase}_1.fastq.gz"
    [[ ! -e $f1 ]] && echo "$f1 not found, exiting" && exit 1
    
    f1out="${outdir}/${dbase}_1.trimmed.fastq.gz"

    BNAME=$outdir/${dbase}.qsub.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; cutadapt -o $f1out $f1 -m 3 -a 'IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT' -a 'IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG' -a  'IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5' -a 'IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -a 'Aseq=TGGCACCCGAGAATTCCA'  -a 'illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT'"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${dbase} --wrap "$cmd"
done


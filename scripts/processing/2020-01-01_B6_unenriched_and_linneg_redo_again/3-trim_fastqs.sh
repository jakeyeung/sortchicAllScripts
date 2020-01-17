#!/bin/sh
# Jake Yeung
# 3-trim_fastqs.sh
#  
# 2019-09-28

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

# WRAP UP
while [[ `qstat | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='4G'
jtime='12:00:00'

for indir in `ls -d $inmain/*BM*`; do
    echo $indir
    bname=$(basename $indir)
    BNAME=$indir/${bname}.qsub.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    f1=$indir/demultiplexedR1.fastq.gz
    f2=$indir/demultiplexedR2.fastq.gz
    f1out=$indir/trimmed.R1.fastq.gz
    f2out=$indir/trimmed.R2.fastq.gz
    [[ -e $f1out ]] && echo "$f1out  found, continuing" && continue
    [[ -e $f2out ]] && echo "$f2out  found, continuing" && continue
    echo "Trimming $f1 $f2"
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; cutadapt -o $f1out -p $f2out $f1 $f2 -m 3 -a 'IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT' -a 'IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG'  -A 'IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5' -a  'IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5' -a 'IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'Aseq=TGGCACCCGAGAATTCCA' -a 'Aseq=TGGCACCCGAGAATTCCA'  -a 'illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT'" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
done

#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --mem=100G
#SBATCH -n 8 -N 1

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/prefetch_outputs"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/prefetch_outputs/fastqs"

cp /hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/fastqs/"$filename"*.fastq.gz .
zcat "$filename"_1.fastq.gz | sed -n '2~4p' > "$filename".index
zcat "$filename"_3.fastq.gz | sed -n '2~4p' > "$filename".cellbc
# rm "$filename"_1.fastq.gz
# rm "$filename"_3.fastq.gz
zcat "$filename"_2.fastq.gz | sed -n '1~4p' | awk -F" " '{print $1}' > "$filename".header
paste -d"." "$filename".header "$filename".index "$filename".cellbc > "$filename".newheader
gunzip "$filename"_2.fastq.gz
gunzip "$filename"_4.fastq.gz
paste -d '\n' <"$filename"_2.fastq "$filename".newheader - - - -| sed '2~5d' | sed '3~4d' | sed '3~3 i +' > "$filename"_2_cbc.fastq 
paste -d '\n' <"$filename"_4.fastq "$filename".newheader - - - -| sed '2~5d' | sed '3~4d' | sed '3~3 i +' > "$filename"_4_cbc.fastq
# rm "$filename"_2.fastq
# rm "$filename"_4.fastq
gzip "$filename"_2_cbc.fastq 
gzip "$filename"_4_cbc.fastq
# rm "$filename".index
# rm "$filename".cellbc
# rm "$filename".header
# rm "$filename".newheader

cutadapt -o "$filename"_2_cbc_trimmed.fastq.gz -p "$filename"_4_cbc_trimmed.fastq.gz "$filename"_2_cbc.fastq.gz "$filename"_4_cbc.fastq.gz -m 3 -a 'IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT' -a 'IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG'  -A 'IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5' -a  'IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5' -a 'IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'Aseq=TGGCACCCGAGAATTCCA' -a 'Aseq=TGGCACCCGAGAATTCCA'  -a 'illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT'

rm "$filename"_2_cbc.fastq.gz
rm "$filename"_4_cbc.fastq.gz

/hpc/hub_oudenaarden/bin/software/bwa-0.7.10/bwa mem -t 8 /hpc/hub_oudenaarden/group_references/ensembl/97/homo_sapiens/primary_assembly_NOMASK_ERCC92_WithLambdaPhage.fa "$filename"_2_cbc_trimmed.fastq.gz "$filename"_4_cbc_trimmed.fastq.gz | samtools view -Sb - > "$filename".bam

rm "$filename"_2_cbc_trimmed.fastq.gz
rm "$filename"_4_cbc_trimmed.fastq.gz

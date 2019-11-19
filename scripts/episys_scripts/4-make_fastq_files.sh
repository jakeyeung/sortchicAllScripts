#!/bin/sh
# Jake Yeung
# 4-make_fastq_files.sh
#  
# 2019-06-25

indir="/Users/yeung/data/scchic/for_episys/sorted_bams_filtered"
outdir="/Users/yeung/data/scchic/for_episys/fastq"

for b in `ls -d $indir/*.bam`; do
		bname=$(basename $b)
		bname=${bname%.*}
		outf1=$outdir/${bname}_R1.fastq
		outf2=$outdir/${bname}_R2.fastq
		bedtools bamtofastq -i $b -fq $outf1 -fq2 $outf2 
done

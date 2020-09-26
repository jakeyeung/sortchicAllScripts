#!/bin/sh
# Jake Yeung
# 0-extend_bed_from_TSS.sh
#  
# 2020-08-17

winsize=25000
winsize=5000
winsize2=`echo "$winsize * 2" | bc -l`
inbed="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.bed.gz"
inbedtmp="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.bed"
outbed="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.${winsize2}.again.bed"

zcat $inbed > $inbedtmp
# exit 0
dos2unix $inbedtmp

# echo "zcat $inbed | awk -v OFS="\t" -v FS="\t" -v dist="$winsize" '{$2=$2-dist; $3=$2+dist} 1' > $outbed"
awk -v OFS="\t" -v FS="\t" -v dist="$winsize" '{tss=$2; strand=$4; $2=tss-dist; $3=tss+dist; $4=$5; $5=strand} 1' $inbedtmp > $outbed
# zcat $inbed.tmp | awk -v dist="$winsize" '{tss=$2; strand=$4; $2=tss-dist; $3=tss+dist; $4=$5; $5=strand} 1' > $outbed
# dos2unix $outbed

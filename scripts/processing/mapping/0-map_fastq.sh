#!/bin/sh
# Jake Yeung
# map_fastq.sh
#  
# 2019-03-20

bwabin="/hpc/hub_oudenaarden/bin/software/bwa-0.7.10/bwa" 
indxbase="/hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_reformat_reg_chr.fa"
f1="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m1-H3K9me3-2_AH3VGVBGX9_S1/demultiplexedR1.fastq.gz"
f2="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m1-H3K9me3-2_AH3VGVBGX9_S1/demultiplexedR2.fastq.gz"

[[ ! -e $indxbase ]] && echo "$indxbase not found, exiting" && exit 1
[[ ! -e $bwabin ]] && echo "$bwabin not found, exiting" && exit 1
[[ ! -e $f1 ]] && echo "$f1 not found, exiting" && exit 1
[[ ! -e $f2 ]] && echo "$f2 not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bwa_output"
cd $outdir
jmem='64G'
jtime='36:00:00'
sname="PZ-BM-m1-H3K9me3-2_AH3VGVBGX9_S1"
outf="$outdir/bwaMapped.bam"
BNAME=$outdir/$sname
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
ncores=16


echo "echo $outdir; $bwabin mem -t 16 $indxbase $f1 $f2 | samtools view -Sb - > $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu

# $bwabin mem -t 8 /hpc/hub_oudenaarden/gene_models/mouse_gene_models/mm10_reformat_reg_chr.fa /hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m2-H3K27me3-2_H2GV2BGX9_S18/demultiplexedR1.fastq.gz /hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m2-H3K27me3-2_H2GV2BGX9_S18/demultiplexedR2.fastq.gz

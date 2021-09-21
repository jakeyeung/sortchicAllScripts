#!/bin/sh
# Jake Yeung
# 2-widen_windows_MouseBM_TSS.sh
# bedfiles are 2bp wide, make them bigger
# 2020-07-28

indir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromRNAseq"

# outdir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromRNAseq/bigger_windows_nochr"
outdir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromRNAseq/bigger_windows_withchr"
[[ ! -d $outdir ]] && mkdir $outdir

winsize=5000

ctypes="Erythroblast HSCs Bcell Neutrophil HighExprs LowExprs"

# for b in `ls -d $indir/*.bed`; do
for ctype in $ctypes; do
    b=${indir}/MouseBM_TSS.${ctype}.bsize_2.bed
    [[ ! -e $b ]] && echo "$b not found, exiting" && exit 1
    # bout=$outdir/MouseBM_TSS.${ctype}.bsize_${winsize}.nochr.bed
    bout=$outdir/MouseBM_TSS.${ctype}.bsize_${winsize}.withchr.bed
    # awk -F $'\t' -v s=$winsize 'BEGIN {OFS=FS} {print $1, $2=$2-5000, $3=$3+5000, $4}' $b | sed 's/^chr//g' > $bout
    awk -F $'\t' -v s=$winsize 'BEGIN {OFS=FS} {print $1, $2=$2-5000, $3=$3+5000, $4}' $b  > $bout
done

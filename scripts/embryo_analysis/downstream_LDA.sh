#!/bin/sh
# Jake Yeung
# downstream_LDA.sh
#  
# 2019-03-28

nn=40
mindist=0.4

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/downstream_LDA_allele_specific.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
robj="/hpc/hub_oudenaarden/avo/scChiC/maria/lda_out_meanfilt.mouse_embryo_K36me3_build95_AS_LDA25.Robj"
[[ ! -e $robj ]] && echo "$robj not found, exiting" && exit 1
# outf="/Users/yeung/data/scchic/embryo_analysis/plots/H3K36me3_embryos.pdf"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/embryo"
outf="$outdir/embryo_LDA_out.nn_${nn}.mindist_${mindist}.pdf"

sweepparams="FALSE"
# cd $wd; Rscript $rs "H3K36me3" $robj $outf $nn $mindist $sweepparams

jmem='32G'
jtime='1:00:00'
BNAME="$outdir/embryo_qsub"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo "cd $wd; Rscript $rs "H3K36me3" $robj $outf $nn $mindist $sweepparams" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
# echo "cd $wd; Rscript $rs "H3K36me3" $robj $outf"

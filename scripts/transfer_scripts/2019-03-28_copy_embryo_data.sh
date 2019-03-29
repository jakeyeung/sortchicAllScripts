#!/bin/sh
# Jake Yeung
# 2019-03-28_copy_embryo_data.sh
# Do some downstream for Maria
# 2019-03-28

inf="t2:/hpc/hub_oudenaarden/avo/scChiC/maria/lda_out_meanfilt.mouse_embryo_K36me3_build95_AS_LDA25.Robj"
outdir="/Users/yeung/data/scchic/from_cluster/embryo"

scp $inf $outdir/.

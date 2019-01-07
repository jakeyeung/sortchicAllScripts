#!/bin/sh
# Jake Yeung
# 2019-01-07_copy_LDA_GREAT_multimark_analyses.sh
# Download LDA across multiple parameters and marks and also their corresponding GREAT downstream 
# 2019-01-07
# 
# 

indir1="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.5_1000"
indir2="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.3_1000"

for indir in $indir1 $indir2; do
    echo $indir
    bname=$(basename $indir)
    outdir="/tmp/$bname"
    scp -r $indir $outdir
done

#!/bin/sh
# Jake Yeung
# 3-run.convert_LDA_to_bins.sh
# Choose best 
# 2019-09-30

n=0
maxjobs=8

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/choose_k_LDA_output.R"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_ZFbonemarrow"

for indir in `ls -d $inmain/lda_outputs*`; do
    outdir=$indir
    for inf in `ls -d $indir/*.Robj`; do
        echo $inf
        bname=$(basename $inf)
        bname=${bname%.*}
        outf=$outdir/$bname.OutObjs.rds  # will add K_choose afterwards
        Rscript $rs $inf $outf&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done


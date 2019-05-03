#!/bin/sh
# Jake Yeung
# 5-run.summarize_activities_all.sh
#  
# 2019-05-01

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/summarize_activities_all.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/mara_output"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/permute_activities_summary"

n=0
maxjobs=16
for indir in `ls -d $inmain/seed_row_*`; do
    outname=$(basename $indir)
    outname=$outname.txt
    outf=$outdir/$outname
    Rscript $indir $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

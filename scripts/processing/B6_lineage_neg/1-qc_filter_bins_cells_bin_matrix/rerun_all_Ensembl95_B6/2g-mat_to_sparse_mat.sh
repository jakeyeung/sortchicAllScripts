#!/bin/sh
# Jake Yeung
# 2g-mat_to_sparse_mat.sh
# Diffiuclt to handle csvs, try rds 
# 2019-05-07

n=0
maxjobs=32
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/first_pipeline/lib/mat_to_sparse_mat.R"

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for indir in `ls -d $inmain/B6-13W1-BM-H3K*-merged`; do
    bname=$(basename $indir)
    tagdir=$indir/tagged
    inf="$tagdir/$bname.filtered.bincounts.slidewin.csv.gz"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outf="$tagdir/$bname.filtered.bincounts.slidewin.rds"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    Rscript $rs $inf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

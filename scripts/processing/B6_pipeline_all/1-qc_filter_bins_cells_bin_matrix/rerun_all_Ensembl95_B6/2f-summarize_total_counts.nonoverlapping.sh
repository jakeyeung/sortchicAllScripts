#!/bin/sh
# Jake Yeung
# 2f-summarize_total_counts.sh
# Summarize total counts : nonoverlapping bins 
# 2019-06-11

# Summ

# rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/first_pipeline/lib/summarize_counts.R"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/1-qc_filter_bins_cells_bin_matrix/lib/summarize_counts.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# inf="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6/B6-13W1-BM-H3K27me3-2-merged/tagged/B6-13W1-BM-H3K27me3-2-merged.filtered.bincounts.slidewin.csv.gz"   # test
# outf="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6/B6-13W1-BM-H3K27me3-2-merged/tagged/B6-13W1-BM-H3K27me3-2-merged.filtered.cellsum.slidewin.csv"
# 
# Rscript $rs $inf $outf


n=0
maxjobs=40

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for indir in `ls -d $inmain/B6-13W1-BM-H3K*-merged`; do
    bname=$(basename $indir)
    tagdir=$indir/tagged
    # inf="$tagdir/$bname.filtered.bincounts.slidewin.csv.gz"
    inf="$tagdir/$bname.filtered.bincounts.offset_0.csv"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outf="$tagdir/$bname.filtered.cellsum.offset_0.csv"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    Rscript $rs $inf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

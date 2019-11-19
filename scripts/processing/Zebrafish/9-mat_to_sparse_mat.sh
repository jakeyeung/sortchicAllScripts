#!/bin/sh
# Jake Yeung
# 5-mat_to_sparse_mat.sh
#  
# 2019-06-17

n=0
maxjobs=7
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/1-qc_filter_bins_cells_bin_matrix/lib/mat_to_sparse_mat_new_slidewin_format.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3910_HD_0/countTables"
outmain="$inmain"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

# winsize=100000
# stepsize=20000
winsize=50000
stepsize=10000

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6
Rscript="/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/Rscript"

for inf in `ls -d $inmain/PZ-ChIC*.csv`; do
    echo $inf
    bname=$(basename $inf)
    bname=${bname%.*}
    outf="$outmain/${bname}.rds"
    # [[ -e $outf ]] && echo "$outf found, continuing" && continue
    $Rscript $rs $inf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

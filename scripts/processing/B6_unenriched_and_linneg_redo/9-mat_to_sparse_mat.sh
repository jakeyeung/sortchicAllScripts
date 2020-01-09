#!/bin/sh
# Jake Yeung
# 5-mat_to_sparse_mat.sh
#  
# 2019-06-17

# WRAP UP
while [[ `qstat | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

n=0
maxjobs=4
rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/mat_to_sparse_mat_new_slidewin_format.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6/countTables"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll/countTables"
outmain="$inmain"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

winsize=100000
stepsize=20000

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6

for inf in `ls -d $inmain/*.csv`; do
    echo $inf
    bname=$(basename $inf)
    bname=${bname%.*}
    outf="$outmain/${bname}.rds"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
    # if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    #     # define maxjobs and n using maxjobsn skeleton
    #     wait # wait until all have finished (not optimal, but most times good enough)
    #     echo $n wait
    # fi
done
# wait

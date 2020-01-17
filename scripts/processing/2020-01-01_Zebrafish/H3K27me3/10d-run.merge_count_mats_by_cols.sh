#!/bin/sh
# Jake Yeung
# 5-run.merge_count_mats_by_cols.sh
#  
# 2019-11-08

jmem='8G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/merge_count_mats_by_cols.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3985/countTables_geneTSS"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3985/countTables_geneTSS/merged"
[[ ! -d $outmain ]] && mkdir $outmain

winsizes="20000 50000 100000"
# marks="H3K4me1 H3K4me3 H3K9me3"
# marks="H3K4me3"  # redo H3K4me3 because broken ln for library 1 was bad
# prefixs="ZFWKM ZFWKMCD41plus"
marks="H3K27me3"  # redo H3K4me3 because broken ln for library 1 was bad
prefixs="ZFWKM"

for winsize in $winsizes; do
    for mark in $marks; do
        for prefix in $prefixs; do
            outbase="PZ-ChIC-${prefix}-${mark}.winsize_${winsize}.merged"
            BNAME="$outmain/$outbase.qsub"
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
            infiles="$inmain/PZ-ChIC-${prefix}-${mark}*.winsize_${winsize}.*.RData"
            outfile=$outmain/PZ-ChIC-${prefix}-${mark}.winsize_${winsize}.merged.RData
            [[  -e $outfile ]] && echo "$outfile found, continuing" && continue
            echo "Rscript $rs --infile $infiles --outfile $outfile" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N Merge_${winsize}_${mark}_${prefix}
            # exit 0 
        done
    done
done


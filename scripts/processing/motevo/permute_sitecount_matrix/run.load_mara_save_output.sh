#!/bin/sh
# Jake Yeung
# run.load_mara_save_output.sh
#  
# 2019-04-10

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/load_mara_save_output.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmem='8G'
jtime='0:30:00'

n=0
maxjobs=8

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permute_summary_activities"
[[ ! -d $outmain ]] && mkdir $outmain

for mark in $marks; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/mara_output"
    nohupmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/nohup_summarize_activites"
    [[ ! -d $nohupmain ]] && mkdir $nohupmain
    [[ ! -d $nohupmain ]] && echo "$nohupmain not found, exiting" && exit 1
    for indir in `ls -d $inmain/seed_row_*`; do
        echo $indir
        bname=$(basename $indir)
        seed=$(echo $bname | cut -d"_" -f3)
        outfile=$outmain/activities_long_${mark}_${seed}.txt

        BNAME=$outmain/load_mara_out_${mark}_${seed}
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        outfilezip=${outfile}.gz

        # echo "Rscript $rs $indir $mark $seed $outfile" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
        [[ -e $outfilezip ]] && echo "$outfilezip found, continuing" && continue
        Rscript $rs $indir $mark $seed $outfile&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done

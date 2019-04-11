#$ -S /bin/bash
#$ -l h_rt=0:30:00
#$ -l h_vmem=1G
#$ -o /hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permute_summary_activities/array.out
#$ -e /hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permute_summary_activities/array.err

# Jake Yeung
# run.load_mara_save_output.sh
#  
# 2019-04-10


if [ "$SGE_TASK_ID" = "undefined" ]; then
  echo "Error. task not defined"
  exit 1
else
  scale=$SGE_TASK_ID
fi

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/load_mara_save_output.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmem='8G'
jtime='0:30:00'

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permute_summary_activities"
[[ ! -d $outmain ]] && mkdir $outmain

i=1
for mark in $marks; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/mara_output"
    nohupmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/nohup_summarize_activites"
    [[ ! -d $nohupmain ]] && mkdir $nohupmain
    [[ ! -d $nohupmain ]] && echo "$nohupmain not found, exiting" && exit 1
    for indir in `ls -d $inmain/seed_row_*`; do
        # echo $indir
        bname=$(basename $indir)
        # echo $bname
        seed=$(echo $bname | cut -d"_" -f3)

        dirarray[$i]=$indir
        seedarray[$i]=$seed
        markarray[$i]=$mark 
        ((i++))
    done
done

mark="${markarray[$SGE_TASK_ID]}"
seed="${seedarray[$SGE_TASK_ID]}"
indir="${dirarray[$SGE_TASK_ID]}"
outfile=$outmain/activities_long_${mark}_${seed}.txt
[[ ! -e $outfile.gz ]] && echo "$outfile.gz not found, exiting" && exit 0

Rscript $rs $indir $mark $seed $outfile

#!/bin/sh
# Jake Yeung
# 1-get_scripts_to_rerun.sh
# Parse log file and get names of scripts tha tneed rerunning 
# 2020-02-18

indir="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/4-run_motevo/log"

cd $indir

prefix="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/H3K27me3/motevo_outputs/split/aa/param_files"
outname="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/4-run_motevo/log/scriptnames_to_rerun.txt"
outname2="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/4-run_motevo/log/scriptdirs_to_rerun.txt"
outnamemerge="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/4-run_motevo/log/scriptpaths_to_rerun.txt"

grep "script_file" eqw_jobs_need_rerun.txt  | tr -d '[:blank:]' | cut -d":" -f2 > $outname
grep "sge_o_workdir" eqw_jobs_need_rerun.txt  | tr -d '[:blank:]' | cut -d":" -f2 > $outname2

paste -d"/" $outname2 $outname > $outnamemerge
rm $outname $outname2



# grep "script_file" eqw_jobs_need_rerun.txt  | tr -d '[:blank:]' | cut -d":" -f2 | sed "s/^/${prefix}/" > $outname
# grep "script_file" eqw_jobs_need_rerun.txt  | tr -d '[:blank:]' | cut -d":" -f2 | sed "s/^/\/hpc\/hub_oudenaarden\/jyeung\/data\/scChiC\/tfbs_output_cluster_BM-AllMerged_Peaks_1000\/H3K27me3\/motevo_outputs\/split\/aa\/param_files"/" > $outname
# add prefix
# echo $outname | sed "s/$prefix//"
# sed "s/$prefix\///" $outname
# sed "s/^/$prefix/" $outname > shared2.txt

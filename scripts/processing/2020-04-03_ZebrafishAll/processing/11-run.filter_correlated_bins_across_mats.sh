#!/bin/sh
# Jake Yeung
# 11-run.filter_correlated_bins_across_mats.sh
#  
# 2020-04-08

# # WRAP UP
# while [[ `qstat | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done


jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/processing/filter_correlated_bins_across_mats.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.filtered_by_counts_TAfrac_var.binfilt"
# outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.filtered_by_counts_TAfrac_var.bincorrfilt"
# inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var"
# outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.corrfilt"
inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-14"
outmain=${inmain}.corrfilt
[[ ! -d $outmain ]] && mkdir $outmain

inf1="$inmain/count_mat.H3K4me1.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.rds"
inf2="$inmain/count_mat.H3K4me3.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.rds"
inf3="$inmain/count_mat.H3K27me3.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.rds"
inf4="$inmain/count_mat.H3K9me3.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.rds"

for inf in $inf1 $inf2 $inf3 $inf4; do
    echo $inf
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
done

outf1="$outmain/count_mat.H3K4me1.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.rds"
outf2="$outmain/count_mat.H3K4me3.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.rds"
outf3="$outmain/count_mat.H3K27me3.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.rds"
outf4="$outmain/count_mat.H3K9me3.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.rds"

outbed=${outmain}/"correlated_bins.bed"
outpdf=${outmain}/"correlated_bins.pdf"

infs="$inf1 $inf2 $inf3 $inf4"
outfs="$outf1 $outf2 $outf3 $outf4"
jnames="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
q="0.97"

BNAME="${outmain}/qsub_correlated_bins"

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infiles $infs -outfiles $outfs -names $jnames -quantilethres $q -badbinsoutbed $outbed -pdfout $outpdf"
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infiles $infs -outfiles $outfs -names $jnames -quantilethres $q -badbinsoutbed $outbed -pdfout $outpdf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N filter_corr_bins

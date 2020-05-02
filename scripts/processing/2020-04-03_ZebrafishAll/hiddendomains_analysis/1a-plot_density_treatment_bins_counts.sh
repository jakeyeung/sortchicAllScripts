#!/bin/sh
# Jake Yeung
# 1a-plot_density_treatment_bins_counts.sh
#  
# 2020-04-24

jmem='4G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/plot_density_treatment_bins.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

thres=20
[[ $thres != [0-9]* ]] && echo "Must be integer: $thres" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs"

for indir in `ls -d $inmain/PZ*`; do
    dname=$(basename $indir)
    inf=${indir}/${dname}_treatment_bins.txt
    outf=${indir}/${dname}_density_plot_bin_counts.pdf

    BNAME=${indir}/plotDens_qsub_${dname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf -threshold $thres" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N plotDens_${dname} -m beas -M j.yeung@hubrecht.eu
done

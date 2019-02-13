#!/bin/sh
# Jake Yeung
# 6-run_MARA_on_LDA_output.sh
# Run MARA on LDA output 
# 2019-02-04

runscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/run_mara_batch_promoters.sh"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

jmark="H3K4me1"
jthres="0.99"
jscale=0
jcenter=0
jbyrow=0
jsuffix=""
# jsuffix="_mm9_v1"

jcenterE="TRUE"

# E="/data/shared/scripts_for_mara/run_mara_on_atger/input_data/exprs_matrix.mat"
# E="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/count_mats_peaks_norm/hiddenDomains_cellmin_100-cellmax_500000-binarize_FALSE-BM_H3K4me1.filt_0.99.txt"

E="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/count_mats_peaks_norm/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_${jmark}.filt_${jthres}.center_${jcenterE}.txt"
Ebase=$(basename $E)
Ebase=${Ebase%.*}
[[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

# N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/sitecount_mats/H3K4me1_sitecount_matrix.norm.txt"
# N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/sitecount_mats/H3K4me1_sitecount_matrix.norm.GeneID.txt"
# N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/sitecount_mats/H3K4me1_sitecount_matrix.norm.GeneID.txt"

# N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/sitecount_mats/H3K4me1_sitecount_matrix.withzeros.txt"
N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/sitecount_mats${jsuffix}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${jbyrow}.txt"
[[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

Nbase=$(basename $N)
Nbase=${Nbase%.*}

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_output/H3K4me1_hiddenDomains/binarize_FALSE_0.995_unnorm_sitecount"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_output/$Ebase-$Nbase-$jsuffix"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_output/H3K4me1_hiddenDomains"

[[ -d $outdir ]] && echo "$outdir found, exiting for safety" && exit 1

[[ ! -d $outdir ]] && mkdir -p $outdir

# [[ -d $outdir ]] && echo "$outdir found, exiting to prevent overwrite" && exit 1
[[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

bash $runscript $E $outdir $N

#!/bin/sh
# Jake Yeung
# 6-run_MARA_on_LDA_output_highK.sh
# Run MARA on LDA output: highK
# 2019-03-13

runscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/run_mara_batch_promoters.sh"

[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

# jmark="H3K4me1"
# jmark="H3K4me3"
jmark="H3K27me3"
# jmark="H3K9me3"
jthres="0.99"
jscale=0
jcenter=0
jbyrow=0
jsuffix=""
jsuffixE="highK"

jcenterE="TRUE"
binarize="TRUE"

E="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/count_mats_peaks_norm/hiddenDomains_cellmin_100-cellmax_500000-binarize_${binarize}-BM_${jmark}.filt_${jthres}.center_${jcenterE}_${jsuffixE}.txt"
Ebase=$(basename $E)
Ebase=${Ebase%.*}
[[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/sitecount_mats${jsuffix}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${jbyrow}.txt"
[[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

Nbase=$(basename $N)
Nbase=${Nbase%.*}

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_output/H3K4me1_hiddenDomains/binarize_FALSE_0.995_unnorm_sitecount"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_output/$Ebase-$Nbase-$jsuffix-$jsuffixE"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_output/H3K4me1_hiddenDomains"

[[ -d $outdir ]] && echo "$outdir found, exiting for safety" && exit 1

[[ ! -d $outdir ]] && mkdir -p $outdir

# [[ -d $outdir ]] && echo "$outdir found, exiting to prevent overwrite" && exit 1
[[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

bash $runscript $E $outdir $N

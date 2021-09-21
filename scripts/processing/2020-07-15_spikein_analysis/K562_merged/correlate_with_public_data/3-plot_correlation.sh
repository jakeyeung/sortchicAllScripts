#!/bin/sh
# Jake Yeung
# 3-plot_correlation.sh
#  
# 2020-11-22

jmem='4G'
jtime='1:00:00'

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary/K562_chipseq_vs_chic_comparison.npz"
jmeth="spearman"
jtype="heatmap"
# jtype="scatterplot"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary"
outf="${outdir}/K562_chipseq_vs_chic_comparison.correlation.jmeth_${jmeth}.jtype_${jtype}.pdf"

BNAME=${outdir}/correlation_qsub_${jmeth}_${jtype}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; plotCorrelation --corData $inf --corMethod $jmeth --whatToPlot $jtype --plotFile $outf --skipZeros --removeOutliers" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

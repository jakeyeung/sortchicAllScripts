#!/bin/sh
# Jake Yeung
# 3-plot_correlation.sh
#  
# 2020-11-22

jmem='4G'
jtime='1:00:00'


# jmeth="spearman"
# jtype="heatmap"
jmeths="pearson spearman"
jtypes="heatmap scatterplot"
# jtype="scatterplot"

jdist="1kb"
jsuffix="TSS_TES"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.log2inputs.${jdist}.${jsuffix}/K562_chipseq_vs_chic_comparison.dist_${jdist}.${jsuffix}.npz"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.log2inputs.${jdist}.${jsuffix}"
for jmeth in $jmeths; do
    for jtype in $jtypes; do

        outf="${outdir}/K562_chipseq_vs_chic_comparison.correlation.jmeth_${jmeth}.jtype_${jtype}.${jdist}.${jsuffix}.log2inputs.pdf"
        BNAME=${outdir}/correlation_qsub_${jmeth}_${jtype}.${jdist}.log2inputs
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; plotCorrelation --corData $inf --corMethod $jmeth --whatToPlot $jtype --plotFile $outf --skipZeros --removeOutliers" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

    done
done


#!/bin/sh
# Jake Yeung
# 8-plot_compare.sh
#  
# 2019-06-25

inf="/Users/yeung/data/scchic/for_episys/bigwig_compare_out/compare_out.npz"

# cor="spearman"
cor="pearson"
# jtype="heatmap"
jtype="scatterplot"
dolog1p=TRUE

outf="/Users/yeung/data/scchic/for_episys/bigwig_compare_out/compare_out_${cor}_${jtype}_log1p_${dolog1p}.pdf"

if [ $dolog1p == TRUE ]
then
	plotCorrelation --corData $inf --corMethod $cor --whatToPlot $jtype -o $outf --removeOutliers --skipZeros --log1p
else
	plotCorrelation --corData $inf --corMethod $cor --whatToPlot $jtype -o $outf --removeOutliers --skipZeros
fi

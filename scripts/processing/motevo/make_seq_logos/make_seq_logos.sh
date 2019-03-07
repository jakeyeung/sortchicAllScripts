#!/bin/sh
# Jake Yeung
# make_seq_logos.sh
#  
# 2019-02-10

rs="/Users/yeung/projects/scchic/scripts/processing/motevo/lib/pwm_to_seqlogo.R"

inmain="/Users/yeung/data/databases/mm10_weight_matrices_v2_split"
outmain="/Users/yeung/data/databases/seqlogos_mm10_weight_matrices_v2_split"

for inf in `ls -d $inmain/*.pwm`; do 
	infbase=$(basename $inf)
	infbase=${infbase%.*}
	outf=$outmain/$infbase.jpg
	Rscript $rs $inf $outf
done

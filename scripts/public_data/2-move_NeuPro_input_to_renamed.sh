#!/bin/sh
# Jake Yeung
# 2-move_NeuPro_input_to_renamed.sh
# Move NeuPro to renamed after converting to bigwig 
# 2019-03-21

inf1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/GSE93127_Input_Neu.bw"
inf2="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/GSE93127_Input_Pro.bw"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/renamed"
outf1=$outdir/input_Neu.bw
outf2=$outdir/input_Pro.bw

[[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1
[[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1

ln -s $inf1 $outf1
ln -s $inf2 $outf2

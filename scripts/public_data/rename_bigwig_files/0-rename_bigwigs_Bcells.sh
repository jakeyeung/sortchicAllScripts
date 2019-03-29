#!/bin/sh
# Jake Yeung
# 0-rename_bigwigs.sh
# Rename neutrophil bigwigs 
# 2019-03-21

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs"
outdir=$indir/renamed
[[ ! -d $outdir ]] && mkdir $outdir

for b in `ls -d $indir/*.bw`; do
    bname=$(basename $b)
    bname=${bname%.*}
    # keep mark and celltype only
    celltype=$(echo $bname | cut -d"_" -f2)

   if [ $celltype == "Mat" ]
   then
     mark=$(echo $bname | cut -d"_" -f4)
     celltype="MatBcell"
     nscores=$(echo $bname | sed 's/[^_]//g' | awk '{ print length }')
     if [ "$nscores" = "4" ]
     then
     	# has rep
        rep=$(echo $bname | cut -d"_" -f5)
     elif [ "$nscores" = "3" ]
     then
        rep="rep1"
    else
        echo "notfound"
     fi

   else
     mark=$(echo $bname | cut -d"_" -f3)
     nscores=$(echo $bname | sed 's/[^_]//g' | awk '{ print length }')
     if [ "$nscores" = "3" ]
     then
     	# has rep
        rep=$(echo $bname | cut -d"_" -f4)
     elif [ "$nscores" = "2" ]
     then
        rep="rep1"
     fi
   fi  
    bnew=${mark}_${celltype}_${rep}
    # echo $nscores $rep
    bout=$outdir/$bnew.bw
    # echo "ln -s $b $bout"
    ln -s $b $bout
done

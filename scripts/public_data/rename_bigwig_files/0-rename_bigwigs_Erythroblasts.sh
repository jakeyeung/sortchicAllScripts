#!/bin/sh
# Jake Yeung
# 0-rename_bigwigs.sh
# Rename neutrophil bigwigs 
# 2019-03-21

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_GenomeRes_2014/bigwig_mm10"
outdir=$indir/renamed
[[ ! -d $outdir ]] && mkdir $outdir

for b in `ls -d $indir/*.bw`; do
    # GSM946527_mm9_wgEncodePsuHistoneMegakaryoH3k04me3BE14halfCd1InputRepSignalRep2.bw
    # or
    # GSM946527_mm9_wgEncodePsuHistoneMegakaryoH3k04me3BE14halfCd1InputSig.bw
    bname=$(basename $b)
    bname=${bname%.*}
    # keep mark and ctype only
    suffix=$(echo $bname | cut -d"_" -f3)
    
    # split suffix  by capital letters into underscores
    strsplit=$(echo $suffix | sed 's/\([A-Z]\)/_\1/g' | sed 's/^ *//; s/ *$//')
    
    ctype=$(echo $strsplit | cut -d"_" -f5)
    mark=$(echo $strsplit | cut -d"_" -f6)  # H3k04me3
    # remove 0
    mark=$(echo $mark | sed 's/0//')
    # change k to K
    mark=$(echo $mark | sed 's/k/K/')
    
    repOrSig=$(echo $strsplit | cut -d"_" -f11)  # Rep or Sig

    if [[ $repOrSig == "Rep" ]]
    then
        # lower case rep
        rep=$(echo $strsplit | cut -d"_" -f13 | tr '[:upper:]' '[:lower:]')
    elif [[ $repOrSig == "Sig" ]]
    then
        rep="input"  # rep name as input
    else
        echo "repOrsSig should be Rep or Sig. Found $repOrSig"
        exit 1
    fi  

    bnew=${mark}_${ctype}_${rep}.bw
    # echo $bname
    # echo $bnew
    bout=$outdir/$bnew
    ln -s $b $bout
done

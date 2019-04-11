#!/bin/sh
# Jake Yeung
# 3-collect_npz_objects_into_one_directory.sh
# Put npz objects into one directory... should check whether there are duplicates?? 
# 2019-03-28

# do public comparisons Neu and MatB
dmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/merged_softlinks_log1p"
[[ ! -d $outdir ]] && mkdir $outdir

for d in $(ls -d $dmain/*build95_log1p_bigwigs); do
    for f in $(ls -d $d/*.npz); do
        ln -s $f $outdir/. 
    done
done


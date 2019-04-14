#!/bin/sh
# Jake Yeung
# 3-collect_npz_objects_into_one_directory.sh
# Put npz objects into one directory... should check whether there are duplicates?? 
# 2019-03-28

# do public comparisons Neu and MatB
dmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/merged_softlinks"
# ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/public_data/lib/npz_to_textfile.py"
# [[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

for d in $(ls -d $dmain/*build95); do
    # dbase=$(basename $d)
    # ctype=$(echo $dbase | cut -d"_" -f1)
    for f in $(ls -d $d/*.npz); do
        # echo $f
        ln -s $f $outdir/. 
        # outf=${outdir}.
        # python $ps $f $outf
    done
done

# # do Lara-Astiaso comparisons
# for D in $(ls -d $dmain/Lara-Astiaso_2-2014_build95); do
#     Dbase=$(basename $D)
#     Ctype=$(echo $Dbase | cut -d"_" -f1)
# done

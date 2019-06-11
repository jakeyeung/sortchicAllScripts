#!/bin/sh
# Jake Yeung
# 4-run.npz_to_text.sh
#  
# 2019-03-29

suffix="_build95_B6_stringent_mergedDir"
ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/public_data/lib/npz_to_textfile.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/merged_softlinks2${suffix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/merged_softlinks_textfile${suffix}"
[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=8
for inf in $(ls -d $indir/*.npz); do
    echo $inf
    inbase=$(basename $inf)
    inbase=${inbase%.*}
    outf=$outdir/$inbase.txt
    python $ps $inf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
echo "Done"

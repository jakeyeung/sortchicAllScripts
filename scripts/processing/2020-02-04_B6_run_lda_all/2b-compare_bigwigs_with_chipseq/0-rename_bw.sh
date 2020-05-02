#!/bin/sh
# Jake Yeung
# 0-rename_bw.sh
#  
# 2020-03-17

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_40.bsize_100"
outmain=${inmain}/renamed

for b in `ls -d $inmain/*.bw`; do
    base=$(basename $b)
    mark=$(echo $base | cut -d"-" -f1)
    bstrip=$(echo $base | cut -d"_" -f2)    
    bstrip2=$(echo $bstrip | cut -d"." -f2)

    echo $bstrip2

    [[ $bstrip2 == "" ]] && echo "File is NA, continuing" && continue

    outdir=${outmain}/${mark}
    [[ ! -d $outdir ]] && mkdir $outdir

    outf=$outdir/${bstrip2}.bw

    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1

    ln -s $b $outf
done

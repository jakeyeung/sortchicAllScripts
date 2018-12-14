#!/bin/sh
# Jake Yeung
# 1-run.filter_umi_mapq.sh
# Filter by MAPQ and UMI within a region  
# 2018-12-14

inmain="/hpc/hub_oudenaarden/Peter/data/VAN2979"
outmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods/VAN2979"
pyscript="/home/hub_oudenaarden/jyeung/projects/histone-mods/filter_umi_mapq.py"
[[ ! -d $outmain ]] && mkdir $outmain

# do BM only
for indir in $(ls -d $inmain/BM*); do
    dname=$(basename $indir)
    outdir=$outmain/$dname
    [[ ! -d $outdir ]] && mkdir $outdir
    for f in $(ls -d $indir/*cell*.bam); do
        fname=$(basename $f)
        fname=${fname%%.*}
        outf=$outdir/$fname.filtered.sorted.bam
        tmpf=/tmp/$fname.filtered.bam
        echo $f
        python $pyscript $f $tmpf $outf --logfile $outdir/$fname.log
        exit 0
        ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done

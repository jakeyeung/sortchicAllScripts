#!/bin/sh
# Jake Yeung
# 1-run.filter_umi_mapq.sh
# Filter by MAPQ and UMI within a region  
# 2018-12-14

n=0
maxjobs=12

inmain="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed"
outmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-PZ-Bl6-BM-Linneg"
pyscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/filter_umi_mapq.py"

[[ ! -e $pyscript ]] && echo "$pyscript not found, exiting" && exit 1
[[ ! -d $outmain ]] && mkdir $outmain

# do BM only
for indir in $(ls -d $inmain/PZ-Bl6-BM-Linneg*); do
    dname=$(basename $indir)
    outdir=$outmain/$dname
    [[ ! -d $outdir ]] && mkdir $outdir
    for f in $(ls -d $indir/bwaMapped.bam); do
        fname=$(basename $f)
        fname=${fname%%.*}
        outf=$outdir/$dname.filtered.sorted.bam
        tmpf=$outdir/$dname.tmp.filtered.bam
        dumpf=$outdir/$dname.dump.filtered.bam
        BNAME=$outdir/$dname
        . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $pyscript $f $tmpf $outf --logfile $outdir/$fname.log --dumpfile $dumpf --no_prefix&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done
wait

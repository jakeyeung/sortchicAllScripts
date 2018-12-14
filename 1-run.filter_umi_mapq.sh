#!/bin/sh
# Jake Yeung
# 1-run.filter_umi_mapq.sh
# Filter by MAPQ and UMI within a region  
# 2018-12-14

jmem='10G'
jtime='1:00:00'

inmain="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed"
outmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods"
pyscript="/home/hub_oudenaarden/jyeung/projects/histone-mods/filter_umi_mapq.py"

[[ ! -d $outmain ]] && mkdir $outmain

# do BM only
for indir in $(ls -d $inmain/PZ-BM*H*H*S*); do
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
        echo $f
        echo "python $pyscript $f $tmpf $outf --logfile $outdir/$fname.log" --dumpfile $dumpf | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
        # python $pyscript $f $tmpf $outf --logfile $outdir/$fname.log --dumpfile $dumpf
        ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
    done
done

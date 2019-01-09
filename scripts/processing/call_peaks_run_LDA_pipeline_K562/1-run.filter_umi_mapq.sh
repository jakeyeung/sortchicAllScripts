#!/bin/sh
# Jake Yeung
# 1-run.filter_umi_mapq.sh
# Filter by MAPQ and UMI within a region  
# 2018-12-14

# # activate your python?
# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh
# conda activate py3

jmem='10G'
jtime='1:00:00'

inmain="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_K562"
pyscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/filter_umi_mapq.py"

[[ ! -e $pyscript ]] && echo "$pyscript not found, exiting" && exit 1

[[ ! -d $outmain ]] && mkdir $outmain

thres=30
marks="H3K27me3 H3K4me3"
genome="hs"  # or mm if bone marrow

# do K562 on two marks H3K27me3 and H3K4me3
for mark in $marks; do
    echo $mark
    for indir in $(ls -d $inmain/PZ-K562-$mark-*); do
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
            echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $pyscript $f $tmpf $outf --logfile $outdir/$fname.log --dumpfile $dumpf --mapq_thres $thres --genome $genome" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
            ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1
        done
    done
done

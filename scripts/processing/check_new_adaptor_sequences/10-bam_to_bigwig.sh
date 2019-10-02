#!/bin/sh
# Jake Yeung
# 10-bam_to_bigwig.sh
#  
# 2019-09-05

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/bam_to_bigwig_mm10.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/bigwigs"
[[ ! -d $outdir ]] && mkdir $outdir

bsize=1000
n=0
maxjobs=3

jmem='4G'
jtime='1:00:00'

for indir in `ls -d $inmain/PZ*`; do
    b="${indir}/tagged/bwaMapped.dedup.sorted.bam"
    [[ ! -e $b ]] && echo "$b not found, exiting" && exit 1
    bbase=$(basename $indir)
    bout=$outdir/$bbase.bw
    BNAME=$outdir/$bbase.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "bash $bs $b $bout $bsize" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bbase
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait



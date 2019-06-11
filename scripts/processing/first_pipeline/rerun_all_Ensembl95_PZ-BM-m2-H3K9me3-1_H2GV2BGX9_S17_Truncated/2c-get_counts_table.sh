#!/bin/sh
# Jake Yeung
# 2c-get_counts_table.sh
#  
# 2019-05-06

n=0
maxjobs=20

bin=100000

offsets="0 20000 40000 60000 80000"

ps="/hpc/hub_oudenaarden/bdebarbanson/internalTools/modularDemultiplexer/taggedBamFileToCountTable.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
for indir in `ls -d $inmain/PZ-BM-m*`; do
    dname=$(basename $indir)
    tagdir=$indir/tagged
    bamname=$dname.filtered.sorted.bam
    inbam=$tagdir/$bamname
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    for offset in $offsets; do
        outname=$dname.filtered.bincounts.offset_$offset.csv
        outf=$tagdir/$outname
        [[ -e $outf ]] && echo "$outf found, continuing" && continue
        . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3;  python $ps $inbam -sampleTags SM -joinedFeatureTags chrom,DS -o $outf -offset $offset -bin $bin --showBinEnd&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done
wait


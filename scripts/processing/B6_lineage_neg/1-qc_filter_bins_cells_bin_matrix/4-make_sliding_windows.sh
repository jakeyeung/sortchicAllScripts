#!/bin/sh
# Jake Yeung
# 4-make_sliding_windows.sh
#  
# 2019-06-17

ps="/home/hub_oudenaarden/jyeung/projects/buys_code/SingleCellMultiOmics/singlecellmultiomics/bamProcessing/bamToCountTable.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-PZ-Bl6-BM-Linneg"

binsize=100000
stepsize=20000

n=0
maxjobs=8

for indir in `ls -d $inmain/PZ-Bl6-BM-Linneg-H3*`; do
    dname=$(basename $indir)
    tagdir=$indir/tagged
    bamname=$dname.filtered.sorted.bam
    inbam=$tagdir/$bamname
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    outname=$dname.filtered.bincounts.${binsize}_${stepsize}.csv
    outf=$tagdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3;  python $ps $inbam -sampleTags SM -joinedFeatureTags reference_name -o $outf -bin $binsize -sliding $stepsize -binTag DS&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

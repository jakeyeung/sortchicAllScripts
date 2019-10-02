#!/bin/sh
# Jake Yeung
# 3-count_TA_frequency.sh
#  
# 2019-06-17

n=0
maxjobs=32

ps="/home/hub_oudenaarden/jyeung/projects/buys_code/SingleCellMultiOmics/singlecellmultiomics/bamProcessing/bamToCountTable.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-PZ-Bl6-BM-Linneg"
inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-PZ-Bl6-BM-Linneg"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
for indir in `ls -d $inmain/PZ-Bl6-BM-Linneg-H3*`; do
    dname=$(basename $indir)
    outdir=$outmain/$dname/tagged
    [[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1 
    tagdir=$indir/tagged
    [[ ! -d $tagdir ]] && echo "$tagdir not found, exiting" && exit 1
    bamname="${dname}.filtered.sorted.bam"
    inbam=$tagdir/$bamname
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    outname=$dname.filtered.RZcounts.csv
    outf=$outdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam -sampleTags SM -featureTags RZ -o $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

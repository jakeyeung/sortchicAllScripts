#!/bin/sh
# Jake Yeung
# 2b-count_RZ_frequencies.sh
# Count RZ frequencies  
# 2019-05-06
# 

n=0
maxjobs=32

ps="/hpc/hub_oudenaarden/bdebarbanson/internalTools/modularDemultiplexer/taggedBamFileToCountTable.py"
# inbam="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6/B6-13W1-BM-H3K27me3-1-merged/tagged/B6-13W1-BM-H3K27me3-1-merged.filtered.sorted.bam"
# outf="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6/B6-13W1-BM-H3K27me3-1-merged/tagged/B6-13W1-BM-H3K27me3-1-merged.filtered.RZcounts.csv"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
for indir in `ls -d $inmain/PZ-BM-m*`; do
    dname=$(basename $indir)
    tagdir=$indir/tagged
    bamname=$dname.filtered.sorted.bam
    inbam=$tagdir/$bamname
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    outname=$dname.filtered.RZcounts.csv
    outf=$tagdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3;python $ps $inbam -sampleTags SM -featureTags RZ -o $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

# taggedBamFileToCountTable.py resorted.featureCounts.bam -sampleTags SM -featureTags RZ -o RZtable.csv
# inf="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6/B6-13W1-BM-H3K27me3-1-merged/tagged/B6-13W1-BM-H3K27me3-1-merged.filtered.sorted.bam"



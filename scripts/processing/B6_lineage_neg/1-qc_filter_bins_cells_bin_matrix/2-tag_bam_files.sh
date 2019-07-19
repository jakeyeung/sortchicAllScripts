#!/bin/sh
# Jake Yeung
# 2-tag_bam_files.sh
#  
# 2019-05-06

suffix="PZ-Bl6-BM-Linneg"
inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-${suffix}"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
ps="/hpc/hub_oudenaarden/bdebarbanson/internalTools/universalBamTagger.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

n=0
maxjobs=32

for indir in `ls -d $inmain/${suffix}*`; do
    dname=$(basename $indir)
    outdir=$indir
    inbam=$indir/$dname.filtered.sorted.bam
    cd $indir
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps --chic --ftag -moleculeRadius 0 $inbam&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

# /hpc/hub_oudenaarden/bdebarbanson/internalTools/universalBamTagger.py --chic --ftag -moleculeRadius 5  -o tagged.bam YOURBAMFILE
# 
# it will convert te read names back to proper read names and assign tags which will allow you to separate molecules and cells in IGV.
# 
# (--chic adds CHIC molecule info, --ftag converts the readname, --molecule radius allows reads from the same molecule to be a bit shifted)
# 
# Best,
# 
# Buys

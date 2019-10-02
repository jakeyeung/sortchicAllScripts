#!/bin/sh
# Jake Yeung
# 1c-check_new_sampnames.sh
#  
# 2019-05-10

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/1-qc_filter_bins_cells_bin_matrix/lib/rename_sampname.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"

n=0
maxjobs=16
for d in `ls -d $indir/B6*`; do
    bname=$(basename $d)
    bamin=$d/$bname.filtered.sorted.bam
    bamout=$d/$bname.filtered.uniqname.sorted.bam
    [[ ! -e $bamout ]] && echo "$bamout not found, exiting" && exit 1
    echo $bamout
    samtools view $bamout | cut -d";" -f13 | uniq -c &
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $bamin $bamout&
    # if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    # 	# define maxjobs and n using maxjobsn skeleton
    #     wait # wait until all have finished (not optimal, but most times good enough)
    #     echo $n wait
    # fi
done
wait

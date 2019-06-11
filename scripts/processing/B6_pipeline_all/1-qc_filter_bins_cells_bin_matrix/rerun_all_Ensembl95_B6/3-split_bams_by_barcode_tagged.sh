#!/bin/sh
# Jake Yeung
# 3-split_bams_by_barcode.sh
#  
# 2019-05-08

n=0
maxjobs=16

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/1-qc_filter_bins_cells_bin_matrix/lib/split_bam_by_barcodes_tagged.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

maindir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
[[ ! -d $maindir ]] && echo "$maindir not found, exiting" && exit 1

for indir in $(ls -d $maindir/B6*-merged); do
    # echo $indir
    fbase=$(basename $indir)
    tagdir=$indir/tagged
    inbam=$tagdir/${fbase}.filtered.uniqname.sorted.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    outdir=$tagdir/split_by_bc
    [[ ! -d $outdir ]] && mkdir $outdir
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam $outdir --add_chr_prefix&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

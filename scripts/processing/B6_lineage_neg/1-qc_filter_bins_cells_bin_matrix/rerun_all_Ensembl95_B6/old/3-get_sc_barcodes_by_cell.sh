#!/bin/sh
# Jake Yeung
# 3-get_sc_barcodes_by_cell.sh
# Get barcodes for each cell (8bp ID and 3bp UMI) 
# 2018-12-14

maindir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"

n=0
maxjobs=16
for indir in $(ls -d $maindir/B6-13W1-BM-H3K*-merged); do
    for b in $(ls -d $indir/*.filtered.sorted.bam); do
        echo $b
        outdir=$indir/barcode_counts
        [[  -d $outdir ]] && echo "$outdir found  continuing" && continue
        [[ ! -d $outdir ]] && mkdir $outdir
        samtools view $b |  grep -oh "bc\:[ACGT]*" | awk '{split($0, a, ":"); print a[2]}' | sort | uniq -c | sort -k1,1nr > $outdir/bc_counts.txt&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done
wait

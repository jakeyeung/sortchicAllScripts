#!/bin/sh
# Jake Yeung
# 2-zip_output_files.sh
# Save space by gzip files 
# 2019-03-26

# mark="H3K9me3"
mark="H3K4me3"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_cluster_build95_CorrPeakFilt/${mark}/motevo_outputs"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

dnames="split merged closestbed_multiple_genes"

n=0
maxjobs=3
for dname in $dnames; do
    jdir=$indir/$dname
    [[ ! -d $jdir ]] && echo "$jdir not found, exiting" && exit 1
    gzip -r $jdir
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done


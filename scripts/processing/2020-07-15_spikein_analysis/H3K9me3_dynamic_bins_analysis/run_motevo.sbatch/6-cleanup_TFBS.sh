#!/bin/sh
# Jake Yeung
# 6-cleanup_TFBS.sh
#  
# 2021-01-30

mark="H3K9me3"
# mark="H3K4me3"
# mark="H3K27me3"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/${mark}"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_ZF-AllMerged2_Peaks_1000/${mark}"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM_dynamic_bins_1000/${mark}"
indir="$inmain/motevo_outputs"
indirf1="$inmain/fasta"
indirf2="$inmain/fastasplit"

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $indirf1 ]] && echo "$indirf1 not found, exiting" && exit 1
[[ ! -d $indirf2 ]] && echo "$indirf2 not found, exiting" && exit 1

cd $inmain

# echo "compresing fasta directories..."
tar -zcvf ${indirf1}.tar.gz $(basename $indirf1)
tar -zcvf ${indirf2}.tar.gz $(basename $indirf2)

dnames="split merged closestbed_multiple_genes"

echo "compressing motevo outputs subdirectories"
n=0
maxjobs=3
cd $indir
for dname in $dnames; do
    jdir=$indir/$dname
    [[ ! -d $jdir ]] && echo "$jdir not found, exiting" && exit 1
    zipf=${jdir}.tar.gz
    echo "compressing $jdir to $zipf"
    # gzip -r $jdir
    echo "tar -zcvf $zipf $(basename $jdir)"
    tar -zcvf $zipf $(basename $jdir)
    # tar -zcvf $zipf $(basename jdir)
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
echo "Done compressing"

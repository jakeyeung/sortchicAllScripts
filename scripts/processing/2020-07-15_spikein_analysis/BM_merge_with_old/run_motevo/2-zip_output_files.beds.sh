#!/bin/sh
# Jake Yeung
# 2-zip_output_files.sh
# Save space by gzip files 
# 2020-02-15

# mark="H3K9me3"
# mark="H3K4me1"
mark="H3K27me3"
# mark="H3K4me3"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/${mark}"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks/${mark}"
indir="$inmain/motevo_outputs"

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# cd $inmain

dname="bed"
echo "compressing motevo outputs subdirectories"
n=0
maxjobs=1
# cd $indir/${dname}
# cd $dname

outname="combined_sites_bed"
jdir=$indir/$dname/$outname
[[ ! -d $jdir ]] && echo "$jdir not found, exiting" && exit 1

wd=$indir/$dname
[[ ! -d $wd ]] && echo "$wd not found, exiting" && exit 1
cd $wd

    [[ ! -d $jdir ]] && echo "$jdir not found, exiting" && exit 1
    zipf=${jdir}.combined.sites.tar.gz
    echo "compressing $jdir to $zipf"
    # gzip -r $jdir
    pwd
    tar -zcvf $zipf $(basename $jdir)
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
echo "Done compressing"

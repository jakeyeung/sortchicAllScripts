#!/bin/sh
# Jake Yeung
# 1-split_bam_by_cluster.K4me1.sh
#  
# 2020-10-21

jmem='16G'
jtime='3:00:00'

ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

infannot="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl/cluster_tables_H3K4me1_BM_all_round2.txt"
[[ ! -e $infannot ]] && echo "$infannot not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"

outdir="${indir}/H3K4me1"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links/H3K4me1"

for f in `ls -d $indir/*K4me1*.bam`; do
    bname=$(basename $f)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $f -annotfile $infannot -outdir $outdir -tagid SM -mapq 40 --add_chr_prefix"
    echo $f
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

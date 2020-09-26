#!/bin/sh
# Jake Yeung
# 2-split_bam_by_chromo.sh
#  
# 2020-09-20

jmem='16G'
jtime='2:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_chromosome.py"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
outdir=${inmain}/spikein_genome
[[ ! -d $outdir ]] && mkdir $outdir

jchromo="J02459.1"
for inf in `ls -d $inmain/*.bam`; do
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.${jchromo}.bam
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf -chromo $jchromo"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

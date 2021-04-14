#!/bin/sh
# Jake Yeung
# 0-split_bam_by_chromosome.sh
#  
# 2020-09-20

jmem='16G'
jtime='2:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/VAN5046_BM/spikein_genome_checks/split_bam_by_chromosome.py"
jchromo="J02459.1"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir=${indir}/split_by_spikein_genome
[[ ! -d $outdir ]] && mkdir $outdir
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams/"

for f in `ls -d $indir/*.bam`; do
    bname=$(basename $f)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.$jchromo.bam
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $f $outf -chromo $jchromo"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

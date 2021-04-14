#!/bin/sh
# Jake Yeung
# 8-split_bam_by_cell_cycle.sh
#  
# 2020-08-16


jmem='16G'
jtime='3:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.cluster_tables"
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/bams_split_by_clusters"
[[ ! -d $outmain ]] && mkdir $outmain

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmark="H3K4me1"

mapq=40



for jmark in $jmarks; do
    inbam="${indir}/K562-EtOH-${jmark}-G1-G2.sorted.tagged.bam"
    annotfile="${annotdir}/${jmark}.cell_cycle_cluster_tables.txt"

    outdir="${outmain}/K562-EtOH-${jmark}-G1-G2"
    [[ ! -d $outdir ]] && mkdir $outdir

    BNAME=${outdir}/${jmark}.sbatch_qlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inbam -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix"

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
done


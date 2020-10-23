#!/bin/sh
# Jake Yeung
# 6-samp_read_umi_summary.spikeinonly.sh
#  
# 2020-09-25

jmem='16G'
jtime='2:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/ChIC_TAPS/count_reads_and_UMIs_in_sample.py"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/spikeins_only"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="${inmain}/read_UMI_summary"
[[ ! -d $outdir ]] && mkdir $outdir

mapq=40
for inbam in `ls -d $inmain/*.tagged.bam`; do

    bname=$(basename $inbam)
    bname=${bname%%.*}

    BNAME=${outdir}/${bname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.mapq_${mapq}.csv
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam $outf -mapq $mapq"

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

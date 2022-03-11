#!/bin/sh
# Jake Yeung
# 5-run.count_reads_from_bed_before_and_after_filtering.sh
#  
# 2021-06-06
# 

jmem='16G'
jtime='1:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/count_reads_from_bed_before_and_after_filtering.py"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3"  # .bed.gz files

inmain2="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/Kaya-Okur_vs_hiddendomains_peaks_filtered"

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/count_before_after_filtering_hiddendomains_filt"
[[ ! -d $outmain ]] && mkdir $outmain

for infile in `ls -d $inmain/*.bed.gz`; do
    fname=$(basename $infile)
    fname=${fname%.*}
    infilefilt=${inmain2}/${fname}.HD_filt.bed
    outfile=${outmain}/${fname}.count_before_after.txt
    [[ ! -e $infilefilt ]] && echo "$infilefilt not found, exiting" && exit 1

    BNAME=${outmain}/${fname}.sbatch_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $infile $infilefilt $outfile"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${fname} --wrap "$cmd"
done


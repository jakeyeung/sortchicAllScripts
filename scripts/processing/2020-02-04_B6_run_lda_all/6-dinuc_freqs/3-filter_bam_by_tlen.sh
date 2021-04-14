#!/bin/sh
# Jake Yeung
# 2-run.bam_PE_length.sh
#  
# 2020-07-02

jmem='16G'
jtime='4:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/6-dinuc_freqs/filter_bam_by_tlen.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17"


# tlenmin=147
tlenmin=190
tlenmax=$tlenmin

outdir="${inmain}/bamFilteredByTlen_${tlenmin}"
[[ ! -d $outdir ]] && mkdir $outdir

for inbam in `ls -d $inmain/*.bam`; do 

    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=${outdir}/tlen_sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outbam=${outdir}/${bname}.tlenfilt.bam
    [[ -e $outbam ]] && echo "$outbam found, continuing" && continue

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam $outbam -length_min $tlenmin -length_max $tlenmax"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=tlen_${bname} --wrap "$cmd"
done


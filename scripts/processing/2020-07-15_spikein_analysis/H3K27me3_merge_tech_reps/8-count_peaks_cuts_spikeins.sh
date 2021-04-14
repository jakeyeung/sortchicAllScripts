#!/bin/sh
# Jake Yeung
# 8-count_peaks_cuts_spikeins.sh
#  
# 2020-11-23

jmem='16G'
jtime='12:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/count_peaks_blacklist.faster.clean.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams"

outdir="${inmain}/counts_in_peaks_vs_nonpeaks_vs_blacklist.BM-AllMerged3.peaks_merged.reseq"
[[ ! -d $outdir ]] && mkdir $outdir

jmark="H3K27me3"

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"

echo $jmark
inbeddir="${inbedmain}/hd_merged.${jmark}.FromR.maxcount_40_60_R"
inbedname="merged.${jmark}.cutoff_analysis.merged.nochr.bed"
inbed=${inbeddir}/${inbedname}
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

for inbam in `ls -d $inmain/*.bam`; do
    inbambase=$(basename $inbam)
    inbambase=${inbambase%.*}
    echo $inbambase

    BNAME=${outdir}/${inbambase}.sbatch_out
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    bname=${inbambase}.HiddenDomains.${ctype}_${mlength}_${jmark}
    outf1=$outdir/${bname}.cuts_in_peaks.csv
    [[ -e $outf1 ]] && echo "outf1 $outf1 found, continuing" && continue
    outf2=$outdir/${bname}.cuts_in_genome.csv
    [[ -e $outf2 ]] && echo "outf2 $outf2 found, continuing" && continue

    BNAME1=$outdir/${bname}.cuts_in_peaks.qsub
    DBASE1=$(dirname "${BNAME1}")
    [[ ! -d $DBASE1 ]] && echo "dbase1 $DBASE2 not found, exiting" && exit 1

    BNAME2=$outdir/${bname}.cuts_in_genome.qsub
    DBASE2=$(dirname "${BNAME2}")
    [[ ! -d $DBASE2 ]] && echo "dbase2 $DBASE2 not found, exiting" && exit 1

    cmd1=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam $outf1 --filterXA -minMQ $mapq --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $inbed -mode cuts_in_peak"
    cmd2=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam $outf2 --filterXA -minMQ $mapq --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -mode cuts_in_chromo"
    echo $cmd1
    echo $cmd2
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME1}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts1_${bname} --wrap "$cmd1"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME2}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts2_${bname} --wrap "$cmd2"
done

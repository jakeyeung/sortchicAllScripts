#!/bin/sh
# Jake Yeung
# 4-run.count_peaks_nonpeaks_blacklist.sh
#  
# 2020-11-07

jmem='8G'
jtime='4:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/K562_merged/peak_calling_hidden_domains/count_peaks_blacklist.faster.clean.py"
# ps2="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/K562_merged/peak_calling_hidden_domains/count_peaksORnonpeaks_blacklist.faster.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
[[ ! -d $inmain ]] && echo "inmain $inmain not found, exiting" && exit 1
outdir="${inmain}/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix"
[[ ! -d $outdir ]] && mkdir $outdir

mapq=40
# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.nochr.SpikeIns.bed"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# chromos="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
# chromos="7"

for jmark in $jmarks; do
    inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt_NoChrPrefix/hiddendomains_outputs/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.1000.cutoff"
    inbed="${inbedmain}/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.1000.cutoff_analysis.bed"

    bname="K562_AllMerged_${jmark}.merged.sorted.tagged.bam"
    inbam=${inmain}/${bname}
    [[ ! -e $inbam ]] && echo "inbam $inbam not found, exiting" && exit 1

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
        # exit 0
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME1}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts1_${bname} --wrap "$cmd1"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME2}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts2_${bname} --wrap "$cmd2"
done

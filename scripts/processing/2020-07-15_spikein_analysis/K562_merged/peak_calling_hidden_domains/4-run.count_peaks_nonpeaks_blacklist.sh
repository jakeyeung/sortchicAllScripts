#!/bin/sh
# Jake Yeung
# 4-run.count_peaks_nonpeaks_blacklist.sh
#  
# 2020-11-07

jmem='8G'
jtime='4:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/K562_merged/peak_calling_hidden_domains/count_peaks_nonpeaks_blacklist.py"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="${inmain}/counts_in_peaks_vs_nonpeaks_vs_blacklist"
[[ ! -d $outdir ]] && mkdir $outdir

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

chromos="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
# chromos="7"

for jmark in $jmarks; do
    inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bams_G1filt_split_by_G1filt_NoChrPrefix/hiddendomains_outputs/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.1000.cutoff"
    inbed="${inbedmain}/K562_AllMerged_${jmark}.merged.sorted.tagged.G1filt.sorted.1000.cutoff_analysis.bed"

    bname="K562_AllMerged_${jmark}.merged.sorted.tagged.bam"
    inbam=${inmain}/${bname}
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1

    for chromo in $chromos; do
        echo $chromo

        outf=$outdir/${bname}.cuts_in_peaks_vs_nonpeaks.${chromo}.csv
        [[ -e $outf ]] && echo "$outf found, continuing" && continue

        BNAME=$outdir/${bname}.cuts_in_peaks_vs_nonpeaks.${chromo}.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam $outf -chromo $chromo --filterXA -minMQ $mapq --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels -bedfile $inbed"
        echo $cmd
        # exit 0
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts_${chromo}_${bname} --wrap "$cmd"
        # exit 0
    done

done

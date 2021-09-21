#!/bin/sh
# Jake Yeung
# 4-run.count_peaks_nonpeaks_blacklist.sh
#  
# 2020-11-07

jmem='8G'
jtime='4:00:00'

# ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/K562_merged/peak_calling_hidden_domains/count_peaks_blacklist.faster.clean.py"
ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/count_peaks_blacklist.faster.clean.py"
# ps2="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/K562_merged/peak_calling_hidden_domains/count_peaksORnonpeaks_blacklist.faster.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/split_by_cluster.MAPQ_40.batch2/merged_bam"
[[ ! -d $inmain ]] && echo "inmain $inmain not found, exiting" && exit 1
outdir="${inmain}/counts_in_peaks_vs_nonpeaks_vs_blacklist.BM-AllMerged3.round2only.eryths"
[[ ! -d $outdir ]] && mkdir $outdir

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"

# jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
jmarks="H3K9me3"

# chromos="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
# chromos="7"

inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_R"
# ctype="Eryths"
ctypegrep="Eryth"  # H3K9me3 is Eryth, H3K27me3 and others is Eryths annoying

for jmark in $jmarks; do
    echo $jmark
    
    if [[ $jmark == "H3K9me3" ]]
    then
        ctype="Eryth"
    else
        ctype="Eryths"
    fi

    echo "Mark: $jmark"
    echo "Ctype: $ctype"

    for inbeddir in `ls -d $inbedmain/maxcount*_${jmark}_*${ctypegrep}*.cutoff`; do
        inbedbase=$(basename $inbeddir)  # maxcount_Rx2500.BM_round1_round2_merged_H3K27me3_DCs.2500.cutoff/
        echo $inbedbase
        # ctype=$(echo $inbedbase | cut -d"." -f2 | cut -d"_" -f6)
        # ctype=$(echo $inbedbase | cut -d"." -f2)
        mlength=$(echo $inbedbase | cut -d"." -f3)
        echo "$ctype $mlength $jmark"
        inbedname="BM_round1_round2_merged_${jmark}_${ctype}.${mlength}.cutoff_analysis.bed"
        inbed=${inbeddir}/${inbedname}
        [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
        # for inbam in `ls -d $inmain/BM_round1_round2_merged_${jmark}*.bam`; do
        for inbam in `ls -d $inmain/BM_round2_all.${jmark}*.bam`; do
            inbambase=$(basename $inbam)
            inbambase=${inbambase%.*}
            echo $inbambase
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
            # echo $cmd1
            # echo $cmd2
            # exit 0
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME1}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts1_${bname} --wrap "$cmd1"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME2}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts2_${bname} --wrap "$cmd2"

        done
    done

done

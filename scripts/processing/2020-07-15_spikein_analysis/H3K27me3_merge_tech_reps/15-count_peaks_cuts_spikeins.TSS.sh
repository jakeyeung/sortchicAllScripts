#!/bin/sh
# Jake Yeung
# 4-run.count_peaks_nonpeaks_blacklist.sh
#  
# 2020-11-07

jmem='8G'
jtime='4:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/count_peaks_blacklist.faster.clean.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/split_by_cluster.MAPQ_40.H3K27me3reseq/H3K27me3/merged_by_ctype"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster"
[[ ! -d $inmain ]] && echo "inmain $inmain not found, exiting" && exit 1
outdir="${inmain}/counts_in_peaks_vs_nonpeaks_vs_blacklist.TSS_10000"
[[ ! -d $outdir ]] && mkdir $outdir

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"

jmarks="H3K27me3"

# inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_R"
# inbedmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/hiddendomains_outputs/hiddendomains_outputs_minlength_500.mincount_-10.FromR.maxcount_10_40_60"
# inbeddir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/hiddendomains_outputs/hiddendomains_outputs_minlength_500.mincount_-10.FromR.maxcount_10_40_60/maxcount_40x500.PZ-BM-H3K27me3-Eryths-merged.500.cutoff"
inbeddir="/hpc/hub_oudenaarden/jyeung/data/databases/refseq"
# inbedname="PZ-BM-H3K27me3-Eryths-merged.500.cutoff_analysis.bed"
inbedname="MmRefseqTss.chromorenamed.10000.again.nochromo.sorted.bed"

ctype="TSS"

for jmark in $jmarks; do
    echo $jmark
    for inbam in `ls -d $inmain/PZ-BM-${jmark}*-merged.bam`; do
            inbambase=$(basename $inbam)
            inbambase=${inbambase%.*}
            echo $inbambase
            bname=${inbambase}.HiddenDomains.${ctype}_${mlength}_${jmark}

            inbed=${inbeddir}/${inbedname}
            [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

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
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME1}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts1_${bname} --wrap "$cmd1"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME2}_%j.log --ntasks=1 --nodes=1 --job-name=MakeCounts2_${bname} --wrap "$cmd2"
    done
done

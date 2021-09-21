#!/bin/sh
# Jake Yeung
# 5-compare_bigwig_with_published.sh
#  
# 2020-02-27
# compare K4me1 with published data
# compare bigwigs for known celltype versus all k4me1 at TSS
# 

jmem='16G'
jtime='1:00:00'

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/bigwigs_to_compare"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/bigwigs_to_compare.log2inputs"
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.log2inputs"
[[ ! -d $outbase ]] && mkdir $outbase
inbws=$(ls -d $indir/*bw | tr '\n' ' ')
blist="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.bed"
outprefix="${outbase}/K562_chipseq_vs_chic_comparison"
outf=${outprefix}.npz
outftxt=${outprefix}.txt

BNAME=${outprefix}.qsub.output

infbadchromos="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/badchromos.txt"
badchromos=$(while read p; do echo $p; done < $infbadchromos | tr '\n' ' ')
echo $badchromos
# badchromos="chrMT chrKI270728.1 chrKI270727.1 chrKI270442.1 chrKI270729.1 chrGL000225.1 chrKI270743.1 chrGL000008.2 chrGL000009.2 chrKI270747.1 chrKI270722.1 chrGL000194.1 chrKI270742.1 chrGL000205.2 chrGL000195.1 chrKI270736.1 chrKI270733.1 chrGL000224.1 chrGL000219.1 chrKI270719.1 chrGL000216.2 chrKI270712.1 chrKI270706.1 chrKI270725.1 chrKI270744.1 chrKI270734.1 chrGL000213.1 chrGL000220.1 "

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary bins -b $inbws -o $outf -bl $blist --outRawCounts $outftxt --chromosomesToSkip $badchromos --binSize 10000" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N bwsum_all
# echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary bins -b $inbws -o $outf -bl $blist --outRawCounts $outftxt --chromosomesToSkip $badchromos"


# ref="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/bigwigs_G1filt_split_by_G1filt.for_chipseq_comparison"
# inbase1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_ENCODE_rerun_from_Buys"
# inbase2="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/renamed"
# 
# blist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
# [[ ! -e $blist ]] && echo "$blist not found, exiting" && exit 1
# 
# 
# chicbwdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster"
# ref="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.BM_giladi.keeptop_1000.unsorted.bed"
# 
# mark="H3K4me1"
# markalt="H3K4me1"
# 
# cd $inbase  # so ls doesn't give full name
# ctypes=$(ls -d *.bw | sed 's/.bw//g' | cut -d"_" -f2 | tr '\n' ' ')
# 
# echo $ctypes
# 
# outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.compare_with_chipseq.${mark}"
# [[ ! -d $outbase ]] && mkdir $outbase
# 
# n=0
# maxjobs=2
# for ctype in $ctypes; do
#     echo $ctype
#     infbase=$inbase/${mark}_${ctype}.bw
#     [[ ! -e $infbase ]] && echo "$infbase not found, exiting" && exit 1
#     infcompares=$(ls -d $chicbwdir/*${markalt}*.bw)
#     infcompares=$(echo $infcompares | tr '\n' ' ')
#     outf="$outbase/${mark}_${ctype}_comparison.npz"
#     [[ -e $outf ]] && echo "$outf  found, continuing" && continue
#     outftxt="$outbase/${mark}_${ctype}_comparison.txt"
#     [[ -e $outftxt ]] && echo "$outftxt found, exiting" && exit 1
# 
#     BNAME=$outbase/${mark}_${ctype}.qsub
#     DBASE=$(dirname "${BNAME}")
#     [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
# 
#     echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary BED-file --BED $ref -b $infbase $infcompares -o $outf -bl $blist --outRawCounts $outftxt" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N bwsum_${ctype}
#     # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary BED-file --bed $ref -b $infbase $infcompares -o $outf -bl $blist --outRawCounts $outftxt"
#     # exit 0
# 
#     # if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
#     #         # define maxjobs and n using maxjobsn skeleton
#     #     wait # wait until all have finished (not optimal, but most times good enough)
#     #     echo $n wait
#     # fi
# done
# 
# echo "Done"

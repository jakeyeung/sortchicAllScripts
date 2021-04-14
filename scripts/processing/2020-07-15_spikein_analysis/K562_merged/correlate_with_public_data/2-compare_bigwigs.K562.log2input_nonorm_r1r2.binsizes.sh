#!/bin/sh
# Jake Yeung
# 5-compare_bigwig_with_published.sh
#  
# 2020-02-27
# compare K4me1 with published data
# compare bigwigs for known celltype versus all k4me1 at TSS
# 

jmem='16G'
jtime='2:00:00'

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/bigwigs_to_compare"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/bigwigs_to_compare.log2input_nonorm_r1r2"
outbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/multibigwigsummary.newbams.log2input_nonorm_r1r2.binsizes"
[[ ! -d $outbase ]] && mkdir $outbase
inbws=$(ls -d $indir/*bw | tr '\n' ' ')
blist="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.bed"

bsizeouts="10000 25000 50000"
# bsizeout="50000"
# bsizeout="100000"

for bsizeout in $bsizeouts; do

    outprefix="${outbase}/K562_chipseq_vs_chic_comparison.bsize_${bsizeout}"
    outf=${outprefix}.npz
    outftxt=${outprefix}.txt

    BNAME=${outprefix}.qsub.output


    infbadchromos="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/compare_with_chipseq_K562/badchromos.txt"
    badchromos=$(while read p; do echo $p; done < $infbadchromos | tr '\n' ' ')
    echo $badchromos
    # badchromos="chrMT chrKI270728.1 chrKI270727.1 chrKI270442.1 chrKI270729.1 chrGL000225.1 chrKI270743.1 chrGL000008.2 chrGL000009.2 chrKI270747.1 chrKI270722.1 chrGL000194.1 chrKI270742.1 chrGL000205.2 chrGL000195.1 chrKI270736.1 chrKI270733.1 chrGL000224.1 chrGL000219.1 chrKI270719.1 chrGL000216.2 chrKI270712.1 chrKI270706.1 chrKI270725.1 chrKI270744.1 chrKI270734.1 chrGL000213.1 chrGL000220.1 "

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; multiBigwigSummary bins -b $inbws -o $outf -bl $blist --outRawCounts $outftxt --chromosomesToSkip $badchromos --binSize $bsizeout" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N bwsum_all

    
done


#!/bin/sh
# Jake Yeung
# 6-run.make_count_matrix_from_bams.sh
# Make count matrix from bams using R 
# 2018-12-19

jmem='16G'
jtime='1:00:00'

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
pvalcutoff="0.3"
mindist="1000"

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $workdir
rs="scripts/processing/make_count_matrix_from_bams.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


for jchip in $marks; do
    # peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.blacklistfilt.broadPeak"
    # peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff/BM_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff_peaks.blacklistfilt.broadPeak"

    # get paths to bams
    bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0"

    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1
    [[ ! -d $bmain ]] && echo "$bmain not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_from_hiddenDomains/count_mats.${pvalcutoff}.${mindist}"
    [[ ! -d $outmain ]] && mkdir $outmain
    outf="$outmain/PZ-BM-${jchip}.merged.NoCountThres.Robj"
    BNAME="$outmain/PZ-BM-${jchip}.merged.NoCountThres"
    tmpf="$outmain/JY_${jchip}_bamlist.out"

    [[ -e $tmpf ]] && echo "$tmpf already exists, skipping $jchip" && continue
    for b in $(ls -d $bmain/PZ*$jchip*/*.sorted.bam); do
        echo $b >> $tmpf
    done
    echo "Temp bam files found in $tmpf"

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "cd $workdir; Rscript $rs $tmpf $peakf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    # Rscript $rs $tmpf $peakf $outf
done


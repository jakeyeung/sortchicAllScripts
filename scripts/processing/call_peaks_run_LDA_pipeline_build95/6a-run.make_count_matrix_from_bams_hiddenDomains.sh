#!/bin/sh
# Jake Yeung
# 6-run.make_count_matrix_from_bams_hiddenDomains.sh
# Make count matrix from bams using R 
# 2019-01-08

jmem='16G'
jtime='1:00:00'

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
mindist="1000"

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $workdir
rs="scripts/processing/make_count_matrix_from_bams.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

cell="BM"

suffix="build95"

bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_${suffix}"
# [[ ! -d $bmain ]] && mkdir $bmain
[[ ! -d $bmain ]] && echo "$bmain not found, exiting" && exit 1
for jchip in $marks; do
    peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output_${suffix}/${cell}_${jchip}_merged.${mindist}.cutoff/${cell}_${jchip}_merged.${mindist}.cutoff_analysis.blacklistfilt.bed"

    # get paths to bams
    # bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0"
    # bmain="count_thres-0_build95/"

    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}"
    [[ ! -d $outmain ]] && mkdir $outmain
    outf="$outmain/PZ-BM-${jchip}.merged.NoCountThres.hiddenDomains.Robj"
    BNAME="$outmain/PZ-BM-${jchip}.merged.NoCountThres.hiddenDomains"
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


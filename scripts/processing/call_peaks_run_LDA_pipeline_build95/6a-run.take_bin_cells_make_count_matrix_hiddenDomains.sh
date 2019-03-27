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

bamnamesdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_peak_analysis_build95"
[[ ! -e $bamnamesdir ]] && echo "$bamnamesdir not found, exiting" && exit 1

jprefix="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_build95/"  # append to beginning of each bamname
[[ ! -d $jprefix ]] && echo "$jprefix not found, exiting" && exit 1

# add forward slash to every back slash
jprefix=$(echo $jprefix | sed 's/\//\\\//g')
echo "Will add prefix to every line:"
echo $jprefix

for jchip in $marks; do
    peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output_${suffix}/${cell}_${jchip}_merged.${mindist}.cutoff/${cell}_${jchip}_merged.${mindist}.cutoff_analysis.blacklistfilt.bed"
    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1

    bamnamesfile=$bamnamesdir/"JY_${jchip}_bamnames.out"
    [[ ! -e $bamnamesfile ]] && echo "$bamnamesfile not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}.cells_from_bin_analysis"
    [[ ! -d $outmain ]] && mkdir $outmain
    outf="$outmain/PZ-BM-${jchip}.merged.NoCountThres.hiddenDomains.Robj"
    BNAME="$outmain/PZ-BM-${jchip}.merged.NoCountThres.hiddenDomains"
    tmpf="$outmain/JY_${jchip}_bamlist.out"

    [[ -e $tmpf ]] && echo "$tmpf already found, exiting" && exit 1

    echo "sed -e s/^/${jprefix}/ $bamnamesfile > $tmpf"
    sed -e "s/^/${jprefix}/" $bamnamesfile > $tmpf

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "cd $workdir; Rscript $rs $tmpf $peakf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    # echo "cd $workdir; Rscript $rs $tmpf $peakf $outf"
done


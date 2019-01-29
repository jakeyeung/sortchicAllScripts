#!/bin/sh
# Jake Yeung
# 4-run.make_count_matrix_from_bams.sh
# Make count matrix from bams using R 
# 2019-01-08

jmem='16G'
jtime='1:00:00'

cell="K562"
marks="H3K27me3 H3K4me1"
mindist="1000"

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $workdir
rs="scripts/processing/make_count_matrix_from_bams.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


for jchip in $marks; do
    peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/hiddenDomains_output_${cell}_round2/${cell}_${jchip}_merged.${mindist}.cutoff/${cell}_${jchip}_merged.${mindist}.cutoff_analysis.blacklistfilt.bed"

    # get paths to bams
    bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc_${cell}_round2/count_thres-0"

    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1
    [[ ! -d $bmain ]] && echo "$bmain not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_${cell}/count_mats.hiddenDomains.${mindist}_round2"
    [[ ! -d $outmain ]] && mkdir $outmain
    outf="$outmain/PZ-${cell}-${jchip}.merged.hiddenDomains.NoCountThres.Robj"
    BNAME="$outmain/PZ-${cell}-${jchip}.merged.hiddenDomains.NoCountThres"
    tmpf="$outmain/JY_${jchip}_bamlist.out"

    [[ -e $tmpf ]] && echo "$tmpf already exists, skipping $jchip" && continue
    for b in $(ls -d $bmain/PZ-${cell}-*${jchip}*/*.sorted.bam); do
        echo $b >> $tmpf
    done
    echo "Temp bam files found in $tmpf"

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "cd $workdir; Rscript $rs $tmpf $peakf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    # Rscript $rs $tmpf $peakf $outf
done


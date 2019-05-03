#!/bin/sh
# Jake Yeung
# 6a-run.take_bin_cells_make_count_matrix_hiddenDomains.sh
# You want to take bin cells because then you have same cells in both activity and bin analysis
# Make count matrix from bams using R 
# 2019-04-15

jmem='16G'
jtime='1:00:00'

n=0
maxjobs=4

# markref="H3K4me1"  # take peak file from one mark, but create count matrix for 4 marks so we can remove bad peaks

# marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
# marks="H3K4me1"
# marks="H3K4me3"
# marks="H3K27me3"
 marks="H3K9me3"
mindist="1000"

# workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"  #  more than 7 gb in tmp files so move to data
workdir="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic"  # more room for Rsubread tmp files
cd $workdir
rs="scripts/processing/make_count_matrix_from_bams.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

cell="BM"

suffix="build95.withchr"
suffix2="CorrPeakFilt"

bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_${suffix}"
# [[ ! -d $bmain ]] && mkdir $bmain
[[ ! -d $bmain ]] && echo "$bmain not found, exiting" && exit 1

bamnamesdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_peak_analysis_${suffix}"
[[ ! -e $bamnamesdir ]] && echo "$bamnamesdir not found, exiting" && exit 1

jprefix="${bmain}/"
# jprefix="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_${suffix}.withchr/"  # append to beginning of each bamname
[[ ! -d $jprefix ]] && echo "$jprefix not found, exiting" && exit 1

# add forward slash to every back slash
jprefix=$(echo $jprefix | sed 's/\//\\\//g')
echo "Will add prefix to every line:"
echo $jprefix

for jmark in $marks; do
    # peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output_${suffix}/${cell}_${jmark}_merged.${mindist}.cutoff/${cell}_${jmark}_merged.${mindist}.cutoff_analysis.blacklistfilt.bed"
    peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_build95/merged_across_clusters_${jmark}/merged_${jmark}.1000.cutoff_analysis.blacklistfilt.CorrPeakFilt.bed"
    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1

    bamnamesfile=$bamnamesdir/"JY_${jmark}_bamnames.out"
    [[ ! -e $bamnamesfile ]] && echo "$bamnamesfile not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}.cells_from_bin_analysis/CorrPeakFilt"
    [[ ! -d $outmain ]] && mkdir $outmain
    outf="$outmain/PZ-BM-${jmark}.merged.NoCountThres.hiddenDomains.${suffix2}.Robj"
    [[ -e $outf ]] && echo "$outf already found, continuing" && continue
    BNAME="$outmain/PZ-BM-${jmark}.merged.NoCountThres.hiddenDomains"
    tmpf="$outmain/JY_${jmark}_bamlist.out"

    if [ ! -e $tmpf ]
    then
        echo "Doing sed"
        sed -e s/^/${jprefix}/ $bamnamesfile > $tmpf
        # sed -e "s/^/${jprefix}/" $bamnamesfile > $tmpf
    else
        [[ -e $tmpf ]] && echo "$tmpf already found, not doing sed"
    fi

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo "cd $workdir; Rscript $rs $tmpf $peakf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    cd $workdir; Rscript $rs $tmpf $peakf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

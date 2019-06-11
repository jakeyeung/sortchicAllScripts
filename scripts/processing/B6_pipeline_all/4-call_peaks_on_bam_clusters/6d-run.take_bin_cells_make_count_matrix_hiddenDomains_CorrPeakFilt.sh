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

# marks="H3K4me3"
marks="H3K27me3 H3K9me3"
mindist="1000"

workdir="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic"  # more room for Rsubread tmp files
cd $workdir
# rs="scripts/processing/make_count_matrix_from_bams.R"
rs="scripts/processing/make_count_matrix_from_bams_DoNotModifySampName.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

suffix="build95_B6"
suffix2="CorrPeakFilt"

bmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6-split_by_bc"
[[ ! -d $bmain ]] && echo "$bmain not found, exiting" && exit 1

bamnamesdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_peak_analysis_${suffix}"
[[ ! -e $bamnamesdir ]] && echo "$bamnamesdir not found, exiting" && exit 1

experiname="B6"
jprefix="${bmain}/"
# jprefix="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_${suffix}.withchr/"  # append to beginning of each bamname
[[ ! -d $jprefix ]] && echo "$jprefix not found, exiting" && exit 1

# add forward slash to every back slash
jprefix=$(echo $jprefix | sed 's/\//\\\//g')
echo "Will add prefix to every line:"
echo $jprefix

for jmark in $marks; do
    # peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output_${suffix}/${cell}_${jmark}_merged.${mindist}.cutoff/${cell}_${jmark}_merged.${mindist}.cutoff_analysis.blacklistfilt.bed"
    peakf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}/merged_across_clusters_${jmark}/merged_${jmark}.1000.cutoff_analysis.blacklistfilt.CorrPeakFilt.bed"
    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1

    bamnamesfile=$bamnamesdir/"JY_${jmark}_bamnames.out"
    [[ ! -e $bamnamesfile ]] && echo "$bamnamesfile not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}.cells_from_bin_analysis/CorrPeakFilt"
    [[ ! -d $outmain ]] && mkdir -p $outmain
    outf="$outmain/${experiname}-${jmark}.merged.NoCountThres.hiddenDomains.${suffix2}.Robj"
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
    cd $workdir; Rscript $rs $tmpf $peakf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

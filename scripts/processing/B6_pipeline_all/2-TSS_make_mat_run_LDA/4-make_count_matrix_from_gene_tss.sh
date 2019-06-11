#!/bin/sh
# Jake Yeung
# 4-make_count_matrix_from_gene_tss.sh
# You want to take bin cells because then you have same cells in both activity and bin analysis
# Make count matrix from bams using R 
# 2019-04-19

jmem='16G'
jtime='1:00:00'

n=0
maxjobs=4

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
mindist="0"
tssdist="20000"

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $workdir
rs="scripts/processing/make_count_matrix_from_bams.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

suffix="build95_B6.withchr"

samp="B6-BM"

# bmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_build95_B6.withchr"
# # [[ ! -d $bmain ]] && mkdir $bmain
# [[ ! -d $bmain ]] && echo "$bmain not found, exiting" && exit 1

# bamnamesdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_peak_analysis_${suffix}"
bamnamesdir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6-bamlist"
[[ ! -e $bamnamesdir ]] && echo "$bamnamesdir not found, exiting" && exit 1

jprefix=""  # no prefix necessary
# jprefix="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_${suffix}.withchr/"  # append to beginning of each bamname
# [[ ! -d $jprefix ]] && echo "$jprefix not found, exiting" && exit 1

# # add forward slash to every back slash
# jprefix=$(echo $jprefix | sed 's/\//\\\//g')
# echo "Will add prefix to every line:"
# echo $jprefix

for jchip in $marks; do
    peakf="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.${tssdist}.bed"
    echo $peakf
    [[ ! -e $peakf ]] && echo "$peakf not found, exiting" && exit 1

    bamnamesfile=$bamnamesdir/"JY_${jchip}_bamlist_goodcellsfilt.out"
    [[ ! -e $bamnamesfile ]] && echo "$bamnamesfile not found, exiting" && exit 1

    # now run Rscript
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.${mindist}_${suffix}_${tssdist}.cells_from_bin_analysis"
    [[ ! -d $outmain ]] && mkdir $outmain
    outf="$outmain/${samp}-${jchip}.merged.NoCountThres.GeneTSS.Robj"
    [[ -e $outf ]] && echo "$outf already found, continuing" && continue
    BNAME="$outmain/${samp}-${jchip}.merged.NoCountThres.GeneTSS"
    tmpf="$outmain/JY_${jchip}_bamlist.out"

    cp $bamnamesfile $tmpf

    # if [ ! -e $tmpf ]
    # then
    #     echo "Doing sed"
    #     sed -e s/^/${jprefix}/ $bamnamesfile > $tmpf
    # else
    #     [[ -e $tmpf ]] && echo "$tmpf already found, not doing sed"
    # fi

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo "cd $workdir; Rscript $rs $tmpf $peakf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err
    cd $workdir; Rscript $rs $tmpf $peakf $outf TRUE&  # TRUE uses pname rather than Chr:Start-End for rownames
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

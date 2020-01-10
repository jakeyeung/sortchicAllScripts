#!/bin/sh
# Jake Yeung
# 6b-filter_counts_blacklist.sh
#  
# 2019-12-25

jmem='48G'
jtime='4:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/filter_good_cells_good_bins.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

dirlh="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/LHcounts"
dircounts="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_otherWinSize"

outmain=${dircounts}/qc_filt
[[ ! -d $outmain ]] && mkdir $outmain
blfile="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.bed.gz"

binsizes="20000 50000"

for binsize in $binsizes; do
    stepsize=$(expr $binsize / 2)

    infcounts1=${dircounts}/H3K4me1-BM_SC-merged.tagged.bsize_${binsize}.step_${stepsize}.countTable.demuxbugfixed.csv
    infcounts2=${dircounts}/H3K4me3-BM_Linneg_SC-merged.tagged.bsize_${binsize}.step_${stepsize}.countTable.demuxbugfixed.csv
    infcounts3=${dircounts}/H3K27me3-BM_Linneg_SC-merged.tagged.bsize_${binsize}.step_${stepsize}.countTable.demuxbugfixed.csv
    infcounts4=${dircounts}/H3K9me3-BM_Linneg_SC-merged.tagged.bsize_${binsize}.step_${stepsize}.countTable.demuxbugfixed.csv

    inf1=${dirlh}/H3K4me1-BM_SC-merged.tagged.LH_counts.csv
    inf2=${dirlh}/H3K4me3-BM_Linneg_SC-merged.tagged.LH_counts.csv
    inf3=${dirlh}/H3K27me3-BM_Linneg_SC-merged.tagged.LH_counts.csv
    inf4=${dirlh}/H3K9me3-BM_Linneg_SC-merged.tagged.LH_counts.csv

    for inf in $inf1 $inf2 $inf3 $inf4 $infcounts1 $infcounts2 $infcounts3 $infcounts4; do
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    done

    outdir=${outmain}/qc_${binsize}_${stepsize}
    [[ ! -d $outdir ]] && mkdir $outdir

    BNAME=${outdir}/qc_${binsize}_${stepsize}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infilerz $inf1 $inf2 $inf3 $inf4 $inf5 -infilecounts $infcounts1 $infcounts2 $infcounts3 $infcounts4 $infcounts5 -names H3K4me1 H3K4me3 H3K27me3 H3K9me3 -outdir $outdir -countcutoff 1000 500 1000 1000 -TAcutoff 0.5 -blfile $blfile --overwrite" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N qc_${binsize}${stepsize}
done

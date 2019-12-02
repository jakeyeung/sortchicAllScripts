#!/bin/sh
# Jake Yeung
# 4-run.project_new_samples_on_LDA_bin.sh
#  
# 2019-09-30

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_stem_cells/1-run_LDA/project_new_samples_on_LDA_bin.R"

jbin="FALSE"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt/projections"
qsubdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt/qsub_outputs"

[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $qsubdir ]] && mkdir $qsubdir

inmat1="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells/PZ-Bl6-BM-StemCells_H3K4me1_matsMergedEnriched_2019-09-29.RData"
inmat2="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells/PZ-Bl6-BM-StemCells_H3K4me1_matsMergedEnrichedNoCellCycle_2019-09-29.RData"

infile="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.OutObjs.RData"

[[ ! -e $infile ]] && echo "$infile not found, exiting" && exit 1
bbase=$(basename $infile)
bbase=${bbase%.*}

jmem='32G'
jtime='24:00:00'

for inmat in $inmat1 $inmat2; do
    inbase=$(echo $bbase | cut -d"_" -f4,5)
    outbase=$(basename $inmat)
    experi=$(echo $outbase | cut -d"_" -f3)

    BNAME=$qsubdir/$inbase.$experi.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo $BNAME

    outfile=$outdir/$inbase.$experi.RData
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $infile $inmat $outfile" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $infile $inmat $outfile --binarizemat"
done


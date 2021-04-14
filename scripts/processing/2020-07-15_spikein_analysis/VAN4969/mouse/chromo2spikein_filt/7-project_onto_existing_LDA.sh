#!/bin/sh
# Jake Yeung
# 7-project_onto_existing_LDA.sh
#  
# 2020-08-22

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

binarize="FALSE"
Kchoose="30"

jmem='32G'
jtime='24:00:00'

inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_H3K4me3_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"

inbase=$(dirname $inlda)
# inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.for_projections/H3K4me3_padded_zeros_for_projections.rds"
inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt/H3K4me3_BM.match_rownames_with_old.rds"

[[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
[[ ! -e $inlda ]] && echo "$inlda not found, continuing" && continue

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_projection_onto_old.chromo2spikeinfilt"
[[ ! -d $outdir ]] && mkdir $outdir

inldabase=$(basename $inlda)
inldabase=${inldabase%.*}

inmatbase=$(basename $inmat)
inmatbase=${inmatbase%.*}

outbase=${inldabase}.x.${inmatbase}

outf="${outdir}/${outbase}.RData"
[[ -e $outf ]] && echo "$outf found, continuing" && continue

basesuffix="${jmark}_${binarize}_${Kchoose}_qsub"
BNAME=${outdir}/${basesuffix}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf"
echo $cmd
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=project_${jmark} --wrap "$cmd"

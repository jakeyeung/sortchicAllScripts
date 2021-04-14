#!/bin/sh
# Jake Yeung
# 2-project_LDA_old_onto_new.sh
#  
# 2020-12-28

jmem='32G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
ldadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_old_new_separate_2020-12-27"
countdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new"

outdir="${ldadir}.projections"
[[ ! -d $outdir ]] && mkdir $outdir

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    inlda="${ldadir}/lda_outputs.count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.new_cells_only.2020-12-27.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.new_cells_only.2020-12-27.K-30.Robj"
    [[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

    inmat="${countdir}/count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.old_cells_only.2020-12-27.rds"
    [[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1

    outbase="BM_project_old_onto_new.${jmark}"
    outf="${outdir}/${outbase}.RData"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=${outdir}/${outbase}.sbatch_output
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=project_${jmark} --wrap "$cmd"
done

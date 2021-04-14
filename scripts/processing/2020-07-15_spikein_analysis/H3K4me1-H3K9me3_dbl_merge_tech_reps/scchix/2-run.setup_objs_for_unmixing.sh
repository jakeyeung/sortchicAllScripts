#!/bin/sh
# Jake Yeung
# 2b-run.setup_objs_for_unmixing.sh
# Visually look for interesting topics, create umap_annotated... RData objects which 
# are input to setting up for unmixing
# 2020-02-05

# WRAP UP
while [[ `squeue -u jyeung | grep MakeClus |  wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/memux_scripts/setup_objs_for_unmixing_topics_general.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


mark1="H3K4me1"
mark2="H3K9me3"
markdbl="H3K4me1xH3K9me3"
suffix="k4_k9_dynamic_regions"
# suffix2="VarFilt_pass2_autosomesOnly.maxcountsfilt"
# experi="mouse_spikein_BMround2all.dbl_common_rows.cellfilt_binfilt"
prefix="ClusterAnnot.lda_and_datmerged.k4_k9_dynamic_regions"
jdate="2021-01-30"

# "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/clustering_outputs.mouse_spikein_BMround2all.dbl_common_rows"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/clustering_outputs.${experi}"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/clustering_outputs.${experi}"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_inputs"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/lda_and_clusters"

# "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/clustering_outputs.mouse_spikein_BMround2all.dbl_common_rows/LouvainAnnot.ldaOut.count_mat.H3K4me1.match_dbl.K-30.RData"

inf1="${inmain}/$prefix.${mark1}.${jdate}.RData"
[[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1

inf2="${inmain}/${prefix}.${mark2}.${jdate}.RData"
[[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1

infdbl="${inmain}/${prefix}.${markdbl}.${jdate}.RData"
[[ ! -e $infdbl ]] && echo "$infdbl not found, exiting" && exit 1

# outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_input/${experi}_${suffix}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/${experi}_${suffix}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_inputs"
[[ ! -d $outdir ]] && mkdir -p $outdir

outprefix="${outdir}/scchix_inputs_clstr_by_celltype_${markdbl}"

BNAME=${outprefix}.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

# echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf1 $inf1 -inf2 $inf2 -infdbl $infdbl -outprefix $outprefix" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N setupobjs_${markdbl}
cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf1 $inf1 -inf2 $inf2 -infdbl $infdbl -outprefix $outprefix"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=SetupObjs_${markdbl} --wrap "$cmd"



#!/bin/sh
# Jake Yeung
# 1-run.1-make_cuts_in_regions.sh
#  
# 2020-12-15

jmem='16G'
jtime='4:00:00'

infbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/bam_dict_dict.pkl"
inftotal="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/total_count_per_cell_dict.pkl"
infnorm="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/normalize_to_counts_dict.pkl"
infgene="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/gene_locations.pkl"
inffeatures="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/features.pkl"

# jgene="Ebf1"
# jgene="Hoxd4"
# jgene="Hoxc9"
# jgene="Zic2"
# jgene="Six3"
jgene="Hoxc6"
percentile="99.5"


ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/plot_cuts_genomic_regions/make_cuts_in_regions.py"

[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
[[ ! -e $inffeatures ]] && echo "$inffeatures not found, exiting" && exit 1
[[ ! -e $infgene ]] && echo "$infgene not found, exiting" && exit 1
[[ ! -e $infnorm ]] && echo "$infnorm not found, exiting" && exit 1
[[ ! -e $inftotal ]] && echo "$inftotal not found, exiting" && exit 1
[[ ! -e $infbam ]] && echo "$infbam not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/scchic_cuts_visualizations_from_sbatch.sigmas"
[[ ! -d $outdir ]] && mkdir $outdir

jradiusleft=20000
# jradiusright=150000
jradiusright=100000
# jmark="H3K4me1"
jmarks="H3K4me1 H3K4me3 H3K27me3"
for jmark in $jmarks; do

    infmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.${jmark}.txt"
    [[ ! -e $infmeta ]] && echo "$infmeta not found, exiting" && exit 1
    BNAME=${outdir}/${jmark}_${jgene}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radiusleft $jradiusleft -radiusright $jradiusright -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile --rserver2hpc_prefix -sigma_cells 0.0001 -sigma 1 -jname $jgene"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radius 100000 -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${jgene} --wrap "$cmd"
    # $cmd

done

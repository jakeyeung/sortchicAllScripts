#!/bin/sh
# Jake Yeung
# 1-run.1-make_cuts_in_regions.sh
#  
# 2020-12-15

jmem='16G'
jtime='4:00:00'

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs_AllFourMarks_again"

infbam="${inmain}/bam_dict_dict.pkl"
inftotal="${inmain}/total_count_per_cell_dict.pkl"
infnorm="${inmain}/normalize_to_counts_dict.pkl"
infgene="${inmain}/gene_locations.pkl"
inffeatures="${inmain}/features.pkl"

# infbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/bam_dict_dict.pkl"
# inftotal="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/total_count_per_cell_dict.pkl"
# infnorm="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/normalize_to_counts_dict.pkl"
# infgene="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/gene_locations.pkl"
# inffeatures="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs2/features.pkl"

# jgene="Hoxd4"
# jgene="Hoxc9"
# jgene="Zic2"
# jgene="Six3"
# jgene="Hoxc6"
# jgene="Hoxc4"
# jgene="Cd34"
# jgene="Erg"
# jgene="Hlf"
# jgene="Kit"
# jgene="Cd34"
# jgene="Hoxa3"
# jgene="Hoxa5"
jgene="Ebf1"
# jgene="S100a8"
# jgene="Cd177"
# jgene="S100a7a"

echo $jgene
percentile="99.5"
# normby="spikein_cuts"
# normby="mean"
# normby="cuts_in_peak"
normby="cuts_total"
jwidth=8
jheight=8
jname="${jgene}_normby_${normby}_jwidth_${jwidth}_height_${jheight}"

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/plot_cuts_genomic_regions/make_cuts_in_regions.py"

[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
[[ ! -e $inffeatures ]] && echo "$inffeatures not found, exiting" && exit 1
[[ ! -e $infgene ]] && echo "$infgene not found, exiting" && exit 1
[[ ! -e $infnorm ]] && echo "$infnorm not found, exiting" && exit 1
[[ ! -e $inftotal ]] && echo "$inftotal not found, exiting" && exit 1
[[ ! -e $infbam ]] && echo "$infbam not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/scchic_cuts_visualizations_from_sbatch.sigmas.FourMarks.ByGene"
[[ ! -d $outdir ]] && mkdir $outdir

outtype="png"
jradiusleft=20000
# jradiusleft=25000
# jradiusright=150000
# jradiusright=25000
jradiusright=60000
# jmark="H3K4me1"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

indirmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
# indirmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"  # for figs 3 and 4

for jmark in $jmarks; do

    # infmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.${jmark}.txt"
    infmeta="${indirmeta}/cell_cluster_table_with_spikeins.${jmark}.2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt"
    # infmeta="${indirmeta}/metadata_batch_corrected.arranged_by_lineage.${jmark}.2021-01-02.txt"
    [[ ! -e $infmeta ]] && echo "$infmeta not found, exiting" && exit 1
    BNAME=${outdir}/${jmark}_${jgene}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radiusleft $jradiusleft -radiusright $jradiusright -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile --rserver2hpc_prefix -sigma_cells 0.0001 -sigma 1 -jname $jname -norm_by $normby -jsep '\t' -colorcname 'clustercol' -outfiletype $outtype -width $jwidth -height $jheight"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radius 100000 -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${jgene} --wrap "$cmd"
    # $cmd

done

#!/bin/sh
# Jake Yeung
# 1-run.1-make_cuts_in_regions.sh
#  
# 2020-12-15

jmem='96'
jtime='4:00:00'

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs_AllFourMarks"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/py_objs_AllFourMarks_again"

infbam="${inmain}/bam_dict_dict.pkl"
inftotal="${inmain}/total_count_per_cell_dict.pkl"
infnorm="${inmain}/normalize_to_counts_dict.pkl"
infgene="${inmain}/gene_locations.pkl"
inffeatures="${inmain}/features.pkl"

# jgene="Ebf1"

# from Peter
# jgene="14:95000000-101000000"
# jgene="8:113000000-119000000"
# Gbe1
# jgene="16:68198090-72420146"

# eryth Hoxc4 loss
# jgene="15:103000000-103050000"

# global HSPC loss
# jgene="5:119450000-119500000"
# jgene="5:119450000-119500000"
# jgene="6:94150000-94200000"

# HSPC-specific
# jgene="16:71792394-73653518"

# Eryth-specific
# jgene="1:157598047-158889362"

# Tead1
# jgene="7:112610950-112951318"

# other global losses
# jgene="1:189149935-189456450"  # Kcnk2
# jgene="9:27328416-27477206"  # Spata19
# jgene="8:44941251-45306635"  # Fat1
# jgene="1:189897214-190816760"  # Prox1
# jgene="17:14219959-14674222"  # Smoc2

# other global losse from Peter
# jgene="4:125000000-132000000"
# jgene="4:126980755-128806842"

# IG regions
# jgene="6:65250725-72603858"  # IGkappa 
# jgene="12:112658112-116530127"  # heavychain

# hox
jgene="6:52100639-52343253"  # hoxa9

# large domain from IGV
# jgene="16:71792394-73653518"
# jgene="1:157598047-158889362"

# normby="spikein_cuts"
# normby="cuts_in_peak"
normby="cuts_total"
# normby="mean"
jgenestr=`echo $jgene | sed 's/\:/_/'`
echo $jgene
echo $jgenestr
# exit 0
percentile="99"
# percentile="99.5"
jwidth=8
jheight=8
smoothing=20

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/plot_cuts_genomic_regions/make_cuts_in_regions.py"

[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
[[ ! -e $inffeatures ]] && echo "$inffeatures not found, exiting" && exit 1
[[ ! -e $infgene ]] && echo "$infgene not found, exiting" && exit 1
[[ ! -e $infnorm ]] && echo "$infnorm not found, exiting" && exit 1
[[ ! -e $inftotal ]] && echo "$inftotal not found, exiting" && exit 1
[[ ! -e $infbam ]] && echo "$infbam not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_jupyter/scchic_cuts_visualizations_from_sbatch.sigmas.fourmarks.shuffled"
[[ ! -d $outdir ]] && mkdir $outdir

jradiusleft=10
# jradiusright=150000
jradiusright=10
# jmark="H3K4me1"
jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmarks="H3K9me3"
# jmarks="H3K9me3"
# jmarks="H3K27me3"
# normby="cuts_in_peak"
# outtype="pdf"

# indirmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
indirmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells"

# jskip=1
# outtype="pdf"

jskip=0
outtype="png"

# sigmacells=0.0001
sigmacells=10
sigmaregion=10

jgenestr="${jgenestr}_normby_${normby}_percentile_${percentile}_width_${jwidth}_height_${jheight}_smooth_${smoothing}_skipheatmap_${jskip}_sigmaregion_${sigmaregion}"

for jmark in $jmarks; do
    # infmeta="${indirmeta}/cell_cluster_table_with_spikeins.${jmark}.2020-12-23.umap_spread.final.txt"
    # infmeta="${indirmeta}/cell_cluster_table_with_spikeins.${jmark}.2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt"
    infmeta="${indirmeta}/metadata_batch_corrected.arranged_by_lineage.shuffled.${jmark}.2021-02-19.txt"
    # infmeta="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/metadata_umap_celltype_cuts.${jmark}.txt"
    [[ ! -e $infmeta ]] && echo "$infmeta not found, exiting" && exit 1
    BNAME=${outdir}/${jmark}_${jgenestr}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    if [ $jskip -eq 1 ]
    then
        echo "Skipping clustermap"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radiusleft $jradiusleft -radiusright $jradiusright -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile --rserver2hpc_prefix -sigma_cells 0.0001 -sigma ${sigmaregion} --is_region -jsep '\t' -colorcname 'clustercol' -outfiletype $outtype -jname $jgenestr -norm_by ${normby} -width $jwidth -height $jheight -trace_sigma ${smoothing} --skip_clustermap"
    else
        echo "Keeping clustermap"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radiusleft $jradiusleft -radiusright $jradiusright -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile --rserver2hpc_prefix -sigma_cells 0.0001 -sigma ${sigmaregion} --is_region -jsep '\t' -colorcname 'clustercol' -outfiletype $outtype -jname $jgenestr -norm_by ${normby} -width $jwidth -height $jheight -trace_sigma ${smoothing}"
    fi

    # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radiusleft $jradiusleft -radiusright $jradiusright -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile --rserver2hpc_prefix -sigma_cells 0.0001 -sigma 1 --is_region -jsep '\t' -colorcname 'clustercol' -outfiletype $outtype -jname $jgenestr -norm_by ${normby} -width $jwidth -height $jheight -trace_sigma ${smoothing} --skip_clustermap"
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py37heatmaps; python $ps -infbam $infbam -inftotalcounts $inftotal -infnorm $infnorm -infgenelocs $infgene -inffeatures $inffeatures -outdir $outdir -radius 100000 -gene $jgene -infmeta $infmeta -mark $jmark -percentile $percentile

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${jgene} --wrap "$cmd"
    # $cmd

done

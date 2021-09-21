#!/bin/sh
# Jake Yeung
# 1c-run.make_clusters_from_topics_mixturemodel.sh
#  
# 2020-04-18

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/dblchic/scripts/processing/1-run_scMEMUX_EtOH_NoTcells_TopBins/make_clusters_from_topics_mixturemodel.R"

jplates="KeepAllPlates"
jmark="K27m3"
# topicskeep="11 25 3 16 21 9"  # look at LDA outputs 
topicskeep="11 25 3 16 9"
jsuffix="TopBins_autosomesOnly"
jsuffixout="${jsuffix}.FewerTopics2"  # try different topics
inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_${jsuffix}"
outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_downstream/clusterbytopics_${jsuffixout}" 

[[ ! -d $outmain ]] && mkdir $outmain
outprefix=${outmain}/clusterbytopics_${jmark}_${jsuffix}_${jplates}
fname="lda_outputs.count_mat.${jmark}.KeepTopBins_500.${jplates}.K-30.binarize.FALSE/ldaOut.count_mat.${jmark}.KeepTopBins_500.${jplates}.K-30.Robj"

inf=${inmain}/${fname}

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

BNAME=$outmain/qsub_${jmark}_${jsuffix}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile ${inf} -outprefix ${outprefix} -topicskeep ${topicskeep} -mark $jmark" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N ${jmark}_${jsuffix} -m beas -M j.yeung@hubrecht.eu


# rs="/home/hub_oudenaarden/jyeung/projects/dblchic/scripts/processing/1-run_scMEMUX_EtOH_NoTcells_TopBins/make_clusters_from_topics_mixturemodel.R"
# 
# jmark="K27m3"
# topicskeep="11 25 3 16 21 9"  # look at LDA outputs 
# jsuffix="TopBins_autosomesOnly"
# inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_${jsuffix}"
# outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_downstream/clusterbytopics_${jsuffix}" 
# [[ ! -d $outmain ]] && mkdir $outmain
# outprefix=${outmain}/clusterbytopics_${jmark}_${jsuffix}
# fname="lda_outputs.count_mat.${jmark}.KeepTopBins_500.PlatesFilt.K-30.binarize.FALSE/ldaOut.count_mat.${jmark}.KeepTopBins_500.PlatesFilt.K-30.Robj"
# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile ${inmain}/${fname} -outprefix ${outprefix} -topicskeep ${topicskeep} -mark $jmark
# 

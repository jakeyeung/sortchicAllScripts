#!/bin/sh
# Jake Yeung
# 1-run.make_clusters_from_LDA_outputs.sh
#  
# 2020-04-04

jmem='16G'
jtime='1:00:00'

inftss="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
[[ ! -e $inftss ]] && echo "$inftss not found, exiting" && exit 1

# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/memux_scripts/make_clusters_from_LDA_outputs.UseClusterColname.R"
rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/memux_scripts/make_clusters_from_csv.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


# jname="BM_EtOH_NoTcells_VarFilt_pass2_autosomesOnly.maxcountsfilt"
jname="mouse_spikein_BMround2all.dbl_common_rows.cellfilt_binfilt"
# inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_${jname}"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_${jname}"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

# outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/clustering_outputs.${jname}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/clustering_outputs.${jname}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/clustering_outputs.${jname}"
[[ ! -d $outdir ]] && mkdir $outdir

# jmarks="K27m3 K9m3 K9m3K27m3"
# # jmarks="K27m3 K9m3 K27m3xK9m3"
# jmarks="K27m3 K9m3 K27m3xK9m3"
jmarks="H3K4me1 H3K9me3"
# jmarks="H3K4me1xH3K9me3"
jsuffix="match_dbl.cellfilt.binfilt"

annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/H3K4me1_H3K9me3_analyses/cluster_tables.withdbl.cellfilt_binfilt"

for jmark in $jmarks; do

  annotf="${annotdir}/cluster_tables_${jmark}_BM_all_round2.txt"
  dname="lda_outputs.count_mat.${jmark}.${jsuffix}.K-30.binarize.FALSE"
  fname="ldaOut.count_mat.${jmark}.${jsuffix}.K-30.Robj"
  inf=${inmain}/${dname}/${fname}
  [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

  bname=${fname%.*}

  outpref=${outdir}/LouvainAnnot.${bname}

  BNAME=$outdir/qsubout.${bname}
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outprefix $outpref -annotfile $annotf -mark $jmark"
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=MakeClusters_${jmark} --wrap "$cmd"
done






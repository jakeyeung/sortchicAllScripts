#!/bin/sh
# Jake Yeung
# 6-run_LDA.sh
#  
# 2020-08-11


jmem='16G'
# jtime='48:00:00'
jtime='48:00:00'


rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables/for_LDA"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/count_tables_from_peaks.from_sitecount_mat/for_LDA"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

# prefix="mouse_spikein_VAN5046"
# prefix="mouse_spikein_VAN5046_varfilt"
prefix="mouse_spikein_BMround2all_old_new_separate_2020-12-27"

outmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins"
[[ ! -d $outmain0 ]] && mkdir $outmain0

outmain="${outmain0}/ldaAnalysisBins_${prefix}"
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

for inf in `ls -d $indir/*.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${bname} --wrap "$cmd"
done

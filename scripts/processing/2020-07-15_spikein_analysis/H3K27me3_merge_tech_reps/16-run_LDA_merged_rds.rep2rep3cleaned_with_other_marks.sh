#!/bin/sh
# Jake Yeung
# 16-run_LDA_merged_rds.rep2rep3cleaned_with_other_marks.sh
#  
# 2021-01-21

jmem='16G'
# jtime='48:00:00'
jtime='60:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/count_tables_from_peaks.from_sitecount_mat/for_LDA"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMfinal.from_bins.10kb"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMfinal.from_bins.10kb"

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

prefix="mouse_spikein_BMround2all_rep2rep3reseq_varfilt_bins_10kb_allmarks"

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

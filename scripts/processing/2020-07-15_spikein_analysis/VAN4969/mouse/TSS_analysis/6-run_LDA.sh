#!/bin/sh
# Jake Yeung
# 6-run_LDA.sh
#  
# 2020-08-11


jmem='32G'
jtime='48:00:00'


rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

prefix="mouse_spikein_TSS"

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

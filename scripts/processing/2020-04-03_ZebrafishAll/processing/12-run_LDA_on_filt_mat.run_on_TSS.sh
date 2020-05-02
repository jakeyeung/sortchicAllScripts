#!/bin/sh
# Jake Yeung
# 1-run_LDA_again_2020-04-03.EtOH_NoTcells.sh
# 2020-04-03

# sleep 600

jmem='16G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"
[[ ! -d $workdir ]] && echo "$workdir not found, exiting" && exit 1

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

# winsize=50000

suffix="TSS"
prefix="ZF_AllMerged.${suffix}"
# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_${winsize}.filtered_by_counts_TAfrac_var"
# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_${winsize}.filtered_by_counts_TAfrac_var.bincorrfilt"
# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.corrfilt"  # 2020-04-12
# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.${suffix}.2020-04-14.corrfilt"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.winsize_50000"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_${prefix}"
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

for inf in `ls -d $indir/*.cellfilt.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    echo $bname

    BNAME=$outmain/$bname.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N LDA.nobin.$bname
done
wait

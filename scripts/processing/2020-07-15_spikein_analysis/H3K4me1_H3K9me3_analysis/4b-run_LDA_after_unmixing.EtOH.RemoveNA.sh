#!/bin/sh
# Jake Yeung
# 2-run_LDA.sh
#  
# 2019-12-06

# # WRAP UP
# while [[ `qstat | grep unmix |  wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

# WRAP UP
while [[ `squeue -u jyeung | grep RunFits | wc -l` > 0 ]]; do
        echo "sleep for 180 seconds"
        sleep 180
done

# echo "Sleep another 60 to be sure"
# sleep 60

# WRAP UP
while [[ `qstat | grep "unmix_ds" | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/dblchic"
[[ ! -d $workdir ]] && echo "$workdir not found, exiting" && exit 1

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

# prefix="EtOH_NoTcells_VarFilt_pass2_autosomesOnly"
# prefix="EtOH_NoTcells_VarFilt_pass2_autosomesOnly.maxcountsfilt"
# suffix="KeepAllPlates"
prefix="mouse_spikein_BMround2all.dbl_common_rows"
suffix="match_dbl"
# markdbl="K27m3-K9m3"
mark1="H3K4me1"
mark2="H3K9me3"
# markdbl="H3K4me1xH3K9me3"

experi=${prefix}_${suffix}
# experilong="MF_BM_EtOH_K27m3-K9m3_split_reads"
# experilong="MF_BM_${prefix}_${suffix}_${markdbl}_split_reads"
experilong="MF_BM_${prefix}_${suffix}_${mark1}-${mark2}_split_reads"

# indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_output_downstream/${experi}/${experilong}"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output_downstream/${prefix}_${suffix}/MF_BM_${prefix}_${suffix}_${mark1}-${mark2}_split_reads"
# matname="MF_${prefix}_${suffix}_clstr_by_louvain.${markdbl}-unmixed_mat.${jmark}.rds"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/unmix_split_${experi}_${experilong}"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/LDA_outputs/unmix_split_${experilong}"
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

for inf in `ls -d $indir/*.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/afterUnmixing.${bname}.binarize.${binarize}"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    echo $bname

    BNAME=$outmain/$bname.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo "cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N LDA.nobin.$bname
    cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDAafter_${bname} --wrap "$cmd"
done
wait

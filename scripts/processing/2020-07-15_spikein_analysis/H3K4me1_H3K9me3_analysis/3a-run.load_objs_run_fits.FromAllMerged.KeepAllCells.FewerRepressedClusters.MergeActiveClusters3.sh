#!/bin/sh
# Jake Yeung
# 2b-run.load_objs_run_fits.from_server.sh
#  
# 2019-08-08

# echo "Sleep an hour"
# sleep 3600

# WRAP UP
while [[ `squeue -u jyeung | grep Setup |  wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# # WRAP UP
# while [[ `qstat | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='32G'
jtime='6:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/load_objs_run_fits.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

workdir="/home/hub_oudenaarden/jyeung/projects/dblchic"
cd $workdir

ncores=24
jmethod="Brent"
# experi="EtOH_NoTcells"
# experi="EtOH_NoTcells_VarFilt_pass2_autosomesOnly"
# experi="EtOH_NoTcells_VarFilt_pass2_autosomesOnly.maxcountsfilt"

# experi="mouse_spikein_BMround2all.dbl_common_rows.cellfilt_binfilt"
# suffix="match_dbl.cellfilt.binfilt"
# experi="mouse_spikein_BMround2all.dbl_common_rows"
# suffix="match_dbl"
experi="SetupObjs_AllMerged"
suffix="UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters3"
removeNA="TRUE"

markdbl="H3K4me1xH3K9me3"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/${experi}_${suffix}_ncores${ncores}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/${experi}_${suffix}_ncores${ncores}.UseClusterColnames"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output/${experi}_${suffix}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_output/${experi}_${suffix}"

# indat="/hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_input/${experi}_${suffix}/${experi}_${suffix}_clstr_by_louvain_${markdbl}.removeNA_FALSE.RData"
# indat="/hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_input/${experi}_${suffix}/${experi}_${suffix}_clstr_by_louvain_${markdbl}.removeNA_FALSE.RData"
# indat="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/${experi}_${suffix}/${experi}_${suffix}_clstr_by_louvain_${markdbl}.removeNA_FALSE.RData"
# indat="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/${experi}_${suffix}/${experi}_${suffix}_clstr_by_louvain_${markdbl}.removeNA_FALSE.RData"
# indat="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/${experi}_${suffix}/${experi}_${suffix}_clstr_by_louvain_${markdbl}.removeNA_${removeNA}.RData"
indat="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_input/${experi}_${suffix}/${experi}_${suffix}.clstr_by_louvain_H3K4me1xH3K9me3.removeNA_${removeNA}.RData"
[[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1

bname=$(basename $indat)
bname=${bname%.*}

[[ ! -d $outdir ]] && mkdir -p $outdir
outdat="${outdir}/unmix_${bname}.RData"

BNAME=$outdir/$bname.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

[[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
[[ -e $outdat ]] && echo "$outdat found, exiting" && exit 1
[[ $ncores != [0-9]* ]] && echo "Must be integer: $ncores" && exit 1

# echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $workdir; Rscript $rs $indat $outdat --ncores $ncores --method $jmethod" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N ${markdbl}_run_fits
cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $workdir; Rscript $rs $indat $outdat --ncores $ncores --method $jmethod"
echo "--cpus-per-task=${ncores}"
echo $cmd
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=${ncores} --job-name=RunFits_${markdbl} --wrap "$cmd"






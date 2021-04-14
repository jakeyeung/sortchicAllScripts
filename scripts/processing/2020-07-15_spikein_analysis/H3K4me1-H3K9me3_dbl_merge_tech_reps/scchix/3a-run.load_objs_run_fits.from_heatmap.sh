#!/bin/sh
# Jake Yeung
# 2b-run.load_objs_run_fits.from_server.sh
#  
# 2019-08-08

# echo "Sleep an hour"
# sleep 3600

# # WRAP UP
# while [[ `squeue -u jyeung | grep Setup |  wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

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

ncores=8
jmethod="Brent"

markdbl="H3K4me1xH3K9me3"
jsuffix2="from_heatmap"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.${jsuffix2}"

indat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_inputs.${jsuffix2}/scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.RData"
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

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $workdir; Rscript $rs $indat $outdat --ncores $ncores --method $jmethod"
echo "--cpus-per-task=${ncores}"
echo $cmd
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=${ncores} --job-name=RunFits_${markdbl} --wrap "$cmd"



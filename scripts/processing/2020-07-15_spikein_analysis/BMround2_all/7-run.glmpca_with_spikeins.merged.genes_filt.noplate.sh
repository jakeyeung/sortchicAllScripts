#!/bin/sh
# Jake Yeung
# 7-run.glmpca_with_spikeins.sh
#  
# 2020-08-11

jmem='16G'
jtime='12:00:00'

prefix="/hpc/hub_oudenaarden"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/VAN4969/glmpca_with_spikeins.R"
outdir="${prefix}/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein_BMround2all"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

inspike="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all/spikein_info_BM_round2_all.txt"
[[ ! -e $inspike ]] && echo "$inspike not found, exiting" && exit 1

indir="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_BM_round2_all"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# jpenalty=10
jpenalty=1
maxiter=5000
minibatch="stochastic"
optimizer="avagrad"
tol="1e-8"
devfilt="5000"

for inf in `ls -d $indir/*.rds`; do
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchout.${jpenalty}.${maxiter}.${minibatch}.${optimizer}.tol_${tol}.devfilt_${devfilt}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.glmpcaout.penalty_${jpenalty}.maxiter_${maxiter}.${minibatch}.${optimizer}.tol_${tol}.devfilt_${devfilt}.RData
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -inspike $inspike -K 30 -penalty $jpenalty -maxIter ${maxiter} -minibatch $minibatch -optimizer $optimizer -tol ${tol} -topndevgenes $devfilt"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=glmpca_${bname} --wrap "$cmd"
done

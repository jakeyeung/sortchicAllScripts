#!/bin/sh
# Jake Yeung
# 7-run.glmpca_with_spikeins.sh
#  
# 2020-08-11

jmem='16G'
jtime='60:00:00'

prefix="/hpc/hub_oudenaarden"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/VAN4969/glmpca_with_spikeins.R"
outspike="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.rds"
outdir="${prefix}/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein.round2.genes_filt"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

indir="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/genes_filt"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

jpenalty=1
maxiter=500
minibatch="stochastic"
optimizer="avagrad"
tol="1e-4"

for inf in `ls -d $indir/K562*.rds`; do
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchout.${jpenalty}.${maxiter}.${minibatch}.${optimizer}.tol_${tol}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.glmpcaout.penalty_${jpenalty}.maxiter_${maxiter}.${minibatch}.${optimizer}.tol_${tol}.RData
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -inspike $outspike -K 30 -penalty $jpenalty -maxIter $maxiter -minibatch $minibatch -optimizer $optimizer -tol $tol"
    # echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=glmpca_${bname} --wrap "$cmd"
    # exit 0
done

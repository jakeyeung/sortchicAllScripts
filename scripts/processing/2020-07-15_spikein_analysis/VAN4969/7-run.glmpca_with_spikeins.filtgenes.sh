#!/bin/sh
# Jake Yeung
# 7-run.glmpca_with_spikeins.sh
#  
# 2020-08-11

jmem='16G'
jtime='24:00:00'

prefix="/hpc/hub_oudenaarden"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/VAN4969/glmpca_with_spikeins.R"
outspike="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData"
outdir="${prefix}/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

indir="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

jpenalty=5
topn=5000

# for inf in `ls -d $indir/*H3K4me3*.rds`; do
for inf in `ls -d $indir/*topn_${topn}.rds`; do
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchout.${jpenalty}.topn_${topn}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.glmpcaout.penalty_${jpenalty}.RData
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -inspike $outspike -K 30 -penalty $jpenalty"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=glmpca_${bname} --wrap "$cmd"
done

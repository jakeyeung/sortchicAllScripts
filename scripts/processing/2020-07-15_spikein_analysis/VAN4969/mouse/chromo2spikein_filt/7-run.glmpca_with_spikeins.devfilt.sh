#!/bin/sh
# Jake Yeung
# 7-run.glmpca_with_spikeins.sh
#  
# 2020-08-11

jmem='16G'
jtime='60:00:00'

prefix="/hpc/hub_oudenaarden"
jsuffix=".chromo2spikeinfilt"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/VAN4969/glmpca_with_spikeins.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inspike="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda${jsuffix}/H3K4me3_BM.spikeins.txt"

outdir="${prefix}/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein${jsuffix}"
[[ ! -d $outdir ]] && mkdir $outdir

jpenalty=5

inf="${prefix}/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda${jsuffix}/H3K4me3_BM.dev_filt.rds"
bname=$(basename $inf)
bname=${bname%.*}

BNAME=${outdir}/${bname}.sbatchout.${jpenalty}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

outf=${outdir}/${bname}.glmpcaout.penalty_${jpenalty}.RData
[[ -e $outf ]] && echo "$outf found, continuing" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -inspike $inspike -K 30 -penalty $jpenalty"
echo $cmd
# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -inspike $inspike -K 30 -penalty $jpenalty
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=glmpca_${bname} --wrap "$cmd"

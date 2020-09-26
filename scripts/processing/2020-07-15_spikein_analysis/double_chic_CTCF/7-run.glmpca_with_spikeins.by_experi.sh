#!/bin/sh
# Jake Yeung
# 7-run.glmpca_with_spikeins.sh
#  
# 2020-08-11

jmem='16G'
jtime='60:00:00'

# prefix="/hpc/hub_oudenaarden"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/double_chic_CTCF/glmpca_with_spikeins.general.R"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA"
inf="${indir}/count_mat_K36me3_CTCF_dbl_TSS.rds"
inspike="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/count_mat_K36me3_CTCF_dbl_CTCFcuts.rds"
outdir="${indir}/glmpcaout"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

jpenalty=1

    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=${outdir}/${bname}.sbatchout.${jpenalty}.by_plate
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/${bname}.glmpcaout.penalty_${jpenalty}.by_plate.RData
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -inspike $inspike -K 30 -penalty $jpenalty"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=glmpca_${bname} --wrap "$cmd"

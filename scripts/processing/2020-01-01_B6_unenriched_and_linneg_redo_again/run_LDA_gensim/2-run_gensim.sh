#!/bin/sh
# Jake Yeung
# 2-run_gensim.sh
#  
# 2019-12-16

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_multicore_gensim.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_B6.2019-12-16/for_gensim"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_python_outputs"
[[ ! -d $outdir ]] && mkdir $outdir

eta="auto"
etashift=0
alpha="auto"
alphashift=0

ntopics=50
ncores=1
npasses=100
nchunks=1152  # 3 plates per chunk

jmem='32G'
jtime='52:00:00'

for inf in `ls -d $inmain/*.mm`; do
    inbase=$(basename $inf)
    inbase=${inbase%%.*}
    outf="$outdir/BM_all.${inbase}.eta_${eta}.alpha_${alpha}.ntopics_${ntopics}.mapq_40.npasses_${npasses}.lda_model"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    BNAME=${outf%.*}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf $ntopics --eta $eta --etashift $etashift --alpha $alpha --alphashift $alphashift --ncores $ncores --npasses $npasses --chunksize $nchunks" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N ${inbase}_${eta}_${alpha}_${ncores}_${ntopics}
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf $ntopics --eta $eta --etashift $etashift --alpha $alpha --alphashift $alphashift --ncores $ncores --npasses $npasses --chunksize $nchunks"
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf $ntopics --eta $eta --etashift $etashift --alpha $alpha --alphashift $alphashift --ncores $ncores --npasses $npasses --chunksize $nchunks"
done

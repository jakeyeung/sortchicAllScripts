#!/bin/sh
# Jake Yeung
# 11-run_LDA_on_rds.sh
#  
# 2019-11-13

jmem='8G'
jtime='48:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_multicore_gensim.py"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

eta="auto"
etashift=0
alpha="auto"
alphashift=0

ntopics=30
ncores=1
npasses=100
nchunks=1152  # 3 plates per chunk

prefix="ZFWKM_allmarks_gensim"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_ZF"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_python_outputs/ldaAnalysisBins_${prefix}"
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1  

for inf in `ls -d $inmain/*.mm`; do
    inbase=$(basename $inf)
    inbase=${inbase%%.*}
    outf="$outdir/ZFWKM_all.${inbase}.eta_${eta}.alpha_${alpha}.ntopics_${ntopics}.mapq_40.npasses_${npasses}.lda_model"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    BNAME=${outf%.*}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf $ntopics --eta $eta --etashift $etashift --alpha $alpha --alphashift $alphashift --ncores $ncores --npasses $npasses --chunksize $nchunks" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N ${inbase}_${eta}_${alpha}_${ncores}_${ntopics}
done


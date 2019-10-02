#!/bin/sh
# Jake Yeung
# run.run_LDA_gensim.sh
# Run LDA on gensim try  
# 2019-07-10

# ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/LDA_python/run_LDA_gensim.py"
ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/LDA_python/run_LDA_multicore_gensim.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/mm_files_for_lda/H3K4me3_BM_sparse_matrix_count_notrans.mm"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/mm_files_for_lda/input_data/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.mm"
ntopics=50
# eta=0.01
eta="auto"
etashift=0
# alpha=0.05  # 50 / ntopics for collapsed gibbs? 
alpha="auto"  # 50 / ntopics
alphashift=0
ncores=1
npasses=250
nchunks=384
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/mm_files_for_lda/more_passes"
outf="${outdir}/H3K4me3_BM_stringent.eta_${eta}.alpha_${alpha}.npasses_${npasses}.ntopics_${ntopics}.nchunks_${nchunks}.ncores_${ncores}.lda_model"

jmem='4G'
jtime='12:00:00'
BNAME=${outf}
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf $ntopics --eta $eta --etashift $etashift --alpha $alpha --alphashift $alphashift --ncores $ncores --npasses $npasses --chunksize $nchunks" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N PZ_BM_stringent_${eta}_${alpha}_${ncores}_${ntopics}

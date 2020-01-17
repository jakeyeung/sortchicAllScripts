#!/bin/sh
# Jake Yeung
# 4-redo_doc_to_topics.sh
# Used old version of LDA gensim the doc to topics didnt include explicit zeros. Fill them in !
# 2019-12-23

jmem='16G'
jtime='2:00:00'
gensimout="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_python_outputs/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16.first_try_doc_to_topics_bug"
gensimin="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_B6.2019-12-16/for_gensim"
ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/write_doc_to_topics_mat_bugfix.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

nohupdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_python_outputs/ldaAnalysisBins_B6BM_All_allmarks.2019-12-16.first_try_doc_to_topics_bug/nohups"
suffix="TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-16.mm"

for f in `ls -d $gensimout/*.lda_model`; do
    fbase=$(basename $f)
    prefix=`echo $fbase | cut -d"." -f2`
    fin=${gensimin}/${prefix}.${suffix}
    [[ ! -e $fin ]] && echo "$fin not found, exiting" && exit 1
    outname="${fbase}.doc_to_topics.zerosbugfixed.csv"
    outf=$gensimout/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$nohupdir/${prefix}.doc2topics.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $fin $f $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
done

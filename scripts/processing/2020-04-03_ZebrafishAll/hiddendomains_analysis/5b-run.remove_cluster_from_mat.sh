#!/bin/sh
# Jake Yeung
# 5b-run.remove_cluster_from_mat.sh
#  
# 2020-08-25

jmem='4G'
jtime='0:15:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/remove_cluster_from_mat.R"

indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.HiddenDomains.imputevarfilt.lessstringent.mapq_40.NewCountFilters/rds_mat_for_LDA"

outdir="${indir}/remove_eryth"
[[ ! -d $outdir ]] && mkdir $outdir

annotdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain"

[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    annot="${annotdir}/cell_to_cluster_table.${jmark}.2020-04-13.txt"
    inbase="${jmark}.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.rds"
    inf="${indir}/${inbase}"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=${inbase%.*}
    outf="${outdir}/${bname}.remove_eryth.rds"

    BNAME=${outdir}/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -annotfile $annot -clstremove eryth -outfile $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -annotfile $annot -clstremove eryth -outfile $outf
done

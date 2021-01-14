#!/bin/sh
# Jake Yeung
# 4b-run_make_exprs_matrix_from_LDA_peaks.sh
#  
# 2020-08-25

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/make_merged_exprs_mat_from_GLMPCA_for_MARA.args.R"
# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/make_merged_exprs_mat_for_MARA.args.R"

jmark="H3K4me1"
# jmark="H3K4me3"
# jmark="H3K27me3"

# jmarks="H3K4me1 H3K4me3 H3K27me3"
# jmarks="H3K4me1 H3K4me3 H3K27me3"
jmarks="H3K4me1 H3K4me3 H3K27me3"
# jmarks="H3K27me3"
# jmarks="H3K4me1 H3K4me3"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun"
 
# jsuffix="bincutoff_0.binskeep_500.byplate.szname_none.niter_500.reorder_rownames.dupfilt"  # use all bins is better

jsuffix="bincutoff_0.binskeep_0.byplate.szname_none.niter_500.reorder_rownames.dupfilt"
# jsuffix="bincutoff_0.binskeep_0.byplate.szname_none.niter_1000.reorder_rownames.dupfilt"

for jmark in $jmarks; do

    echo $jmark
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/${jmark}/mara_input"
    fname="glmpca.${jmark}.${jsuffix}.RData"
    # fname="glmpca.${jmark}.bincutoff_0.binskeep_0.byplate.szname_none.reorder_rownames.dupfilt.RData"
    # fname="glmpca.${jmark}.bincutoff_0.binskeep_0.byplate.szname_none.reorder_rownames.dupfilt.RData"
    # fname="glmpca_${jmark}.RData"

    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks"
    [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

    inf=${indir}/${fname}
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    bname="BM_${jmark}.BM_AllMerged3.glmpca_plate.${jsuffix}"

    outdir="${outmain}/count_mats_peaks_from_GLMPCA_plate"
    [[ ! -d $outdir ]] && mkdir $outdir
    outf="${outdir}"/${bname}.cleanuprows.same_annot_file.txt
    outpdf="${outdir}"/${bname}.cleanuprows.same_annot_file.pdf

    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=${outdir}/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ${BNAME}
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -outpdf $outpdf" 
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -outpdf $outpdf --AddChr" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu

done



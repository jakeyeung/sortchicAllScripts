#!/bin/sh
# Jake Yeung
# 1-run_glmpca_using_LDA_input.sh
#  
# 2020-11-17

jmem='32G'
jtime='24:00:00'

# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init.R"
rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.R"

jmarks="H3K27me3"
# jmarks="H3K27me3"
# jmark="H3K4me1"
# jcname="cell.var.within.sum.norm.log2.CenteredAndScaled"
platename="plate"
# platename="jrep"
szname="none"

# it greps these rownames should be ok, rownames expected in form 1:344-6000;chr1:344-6000
# jchromos="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"

# bincutoff=10000
# bincutoff=5000
bincutoff=0
binskeepvec="0"
binskeep=0
# binskeep=150
# niter=500
nitervec="500 1000"
# nitervec="1000"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep1rep2rep3reseq.peaks.varfilt"
[[ ! -d $outdir ]] && mkdir $outdir

# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2"
annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old"

jsuffix="_peaks"
# jsuffix2="TES"

for binskeep in $binskeepvec; do
    for niter in $nitervec; do
        echo $binskeep
        echo $niter

        for jmark in $jmarks; do
            # infpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat/lda_outputs.count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.K-30.Robj"
            # infpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.K-30.Robj"
            # infpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq${jsuffix}/lda_outputs.PZ-BM-rep3-${jmark}-rep2rep3reseq.${jsuffix2}.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-${jmark}-rep2rep3reseq.${jsuffix2}.K-30.Robj"
            # infpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_${jmark}.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_${jmark}.cleaned.varfilt_2.K-30.Robj"

            # infpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt_peaks/lda_outputs.PZ-BM-rep2rep3reseq-${jmark}.peaks.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep2rep3reseq-${jmark}.peaks.varfilt.K-30.Robj"
            infpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep1rep2rep3reseq_varfilt_2020-12-15/lda_outputs.mat_${jmark}_rep1rep2rep3reseq.peaks.K-30.binarize.FALSE/ldaOut.mat_${jmark}_rep1rep2rep3reseq.peaks.K-30.Robj"

            outbase=${outdir}/glmpca.${jmark}.bincutoff_${bincutoff}.binskeep_${binskeep}.platename_${platename}.szname_${szname}.niter_${niter}.reorder_rownames.dupfilt.suffix${jsuffix}
            # check output doesnt exist
            outcheck="${outbase}.RData"
            [[ -e $outcheck ]] && echo "$outcheck found, continuing" && continue
            BNAME=${outdir}/glmpca.${jmark}.bincutoff_${bincutoff}.binskeep_${binskeep}.platename_${platename}.szname_${szname}.niter_${niter}.reorder_rownames.dupfilt.suffix${jsuffix}
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

            # annotname="cell_cluster_table_with_spikeins.${jmark}.2020-11-18.dupfilt.txt"
            annotname="BM_rep2rep3reseq_${jmark}.cleaned.varfilt_2.metadata.txt"
            annotname="mat_${jmark}_rep1rep2rep3reseq.metadata.txt"

            annotf=${annotdir}/${annotname}

            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infpeaks $infpeaks -infmeta $annotf -outbase $outbase -platecname $platename -sizefactorcname $szname -bincutoff $bincutoff -binskeep $binskeep -niter $niter"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
        done

    done
done



#!/bin/sh
# Jake Yeung
# 1-run_glmpca_using_LDA_input.sh
#  
# 2020-11-17

jmem='32G'
jtime='24:00:00'

# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init.R"
# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.R"
# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.from_project.R"
rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.from_imputed.R"

# jmarks="H3K27me3"
# jmarks="H3K4me3"
jmarks="H3K4me1 H3K4me3 H3K9me3"
# jcname="cell.var.within.sum.norm.log2.CenteredAndScaled"
# platename="plate"
platename="jrep"
szname="none"

# it greps these rownames should be ok, rownames expected in form 1:344-6000;chr1:344-6000
# jchromos="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"

# bincutoff=10000
# bincutoff=5000
# bincutoff=0
# binskeepvec="0 1000"
# binskeep=0
# binskeep=150
# niter=500
nitervec="100"
# nitervec="1000"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/H3K27me3_rep1rep2rep3reseq.peaks.varfilt"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/BM_from_matadj"
[[ ! -d $outdir ]] && mkdir $outdir

# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq/varfilt_2"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old"
annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new"

# indirpeaks="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tm_results_from_projs"
indirimputed="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_batch_correction_output"
indircounts="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_batch_correction_output"

# jsuffix="_old_to_new"
# jsuffix2="TES"

for niter in $nitervec; do
    echo $niter

    for jmark in $jmarks; do
        # fnamepeaks="tm_result${jsuffix}.${jmark}.2020-12-28.RData"
        fnameimputed="mat_adj.${jmark}.from_LDA.binskeep_1000.rds"
        # fnamecounts="mat_adj.${jmark}.from_LDA.binskeep_1000.rds"
        fnamecounts="count_mat.${jmark}.from_LDA.rds"
        infimputed="${indirimputed}/$fnameimputed"
        infcounts="${indircounts}/$fnamecounts"

        outbase=${outdir}/glmpca.${jmark}.from_matadj.platename_${platename}.szname_${szname}.niter_${niter}
        # check output doesnt exist
        outcheck="${outbase}.RData"
        [[ -e $outcheck ]] && echo "$outcheck found, continuing" && continue
        BNAME=${outdir}/glmpca.${jmark}.from_matadj.platename_${platename}.szname_${szname}.niter_${niter}.sbatch
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        annotname="count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.metadata.2020-12-27.txt"
        annotf=${annotdir}/${annotname}
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infimputed $infimputed -infcounts $infcounts -infmeta $annotf -outbase $outbase -platecname $platename -sizefactorcname $szname -niter $niter -jntopics 30"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
    done
done

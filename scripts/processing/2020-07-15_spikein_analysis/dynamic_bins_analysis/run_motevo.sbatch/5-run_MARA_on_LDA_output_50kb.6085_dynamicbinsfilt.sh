#!/bin/sh
# Jake Yeung
# 6-run_MARA_on_LDA_output_highK.sh
# Run MARA on LDA output: highK
# 2019-03-13

# runscript="/home/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/run_mara_batch_promoters.sh"
runscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/run_mara_batch_promoters.sh"
[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

# jmarks="H3K4me3"
# jmarks="H3K4me3"
# jmarks="H3K4me3 H3K27me3"
# jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
jmarks="H3K27me3"
# jmarks="H3K9me3"
# jmarks="H3K4me3"
# jmarksitecounts="H3K9me3"  # bins from H3K9me3 dynamic bins
jmarksitecounts="H3K27me3"  # bins from H3K9me3 dynamic bins

jscale=0
jcenter=0
jbyrow=0

jcenterE="TRUE"

n=0
maxjobs=1

jmem='96G'
jtime='6:00:00'

# jbinskeep=250
# jtopics=30
# pathprefix="/home/jyeung/hpc"
pathprefix="/hpc/hub_oudenaarden/jyeung/data"
# marainputdir="${pathprefix}/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000"
# marainputdir="${pathprefix}/scChiC/mara_analysis_BM-AllMerged3_Peaks"

jsuffix2="dynamic_bins_top_6085"
marainputdir="${pathprefix}/scChiC/mara_analysis_${jsuffix2}"
# Edir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks/count_mats_peaks_k9dynamicbins"
Edir="${marainputdir}/count_mats_50kb_bins"

jmem='96G'
jtime='3:00:00'

for jmark in $jmarks; do
    # Efname="countmat_PZ_fromHiddenDomains_${jmark}.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_${jbinskeep}.ntopics_${jtopics}.2020-02-11.txt"
    # Efname="ldaOut.${jmark}.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.txt"
    # Efname="ldaOut.count_mat_by_peaks.${jmark}.overlap_dynamic_k9_bins.K-30.keepNbins_0.txt"
    # Efname="ldaOut.count_mat_k9_dynamic_bins_50kb.${jmark}.2021-01-28.K-30.keepNbins_0.txt"
    # Efname="ldaOut.count_tables_merged.${jmark}.DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.${jmark}.2021-02-15.txt.K-30.keepNbins_0.txt"
    Efname="ldaOut.count_tables_merged.${jmark}.DE_bins_all_marks_top_6085_dists_to_TSS.annot_table.${jmark}.2021-02-15.txt.K-30.keepNbins_0.withchr2.txt"
    # E="${marainputdir}/${jmark}/mara_input/count_mats_peaks_norm/${Efname}"
    E=${Edir}/${Efname}
    Ebase=$(basename $E)
    Ebase=${Ebase%.*}
    [[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

    # Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.txt"
    # E="${marainputdir}/${jmark}/mara_input/count_mats_peaks_norm/${Efname}"
    E=${Edir}/${Efname}
    Ebase=$(basename $E)
    Ebase=${Ebase%.*}
    [[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

    # Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.txt"
    # Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.txt"
    Nfname="dynamic_bins_50kb__motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmarksitecounts}.YyMerged.tab.txt"
    N="${marainputdir}/${jmarksitecounts}/mara_input/sitecount_mats/${Nfname}"
    [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

    Nbase=$(basename $N)
    Nbase=${Nbase%.*}

    outdir="${marainputdir}/${jmark}/mara_output/${Ebase}-${Nbase}.fromnode"
    [[ -d $outdir ]] && echo "$outdir found, continuing for safety" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

    BNAME=$outdir/${jmark}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "Running MARA for $jmark"
    echo "bash $runscript $E $outdir $N&"
    # bash $runscript $E $outdir $N&
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; bash $runscript $E $outdir $N" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MARA_${jmark}
    # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; bash $runscript $E $outdir $N"
    # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name="${Ebase}x${Nbase}" --wrap "$cmd"
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; bash $runscript $E $outdir $N
    # echo "bash $runscript $E $outdir $N" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done

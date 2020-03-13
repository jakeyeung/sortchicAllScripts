#!/bin/sh
# Jake Yeung
# 6-run_MARA_on_LDA_output_highK.sh
# Run MARA on LDA output: highK
# 2019-03-13

# runscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/run_mara_batch_promoters.sh"
# runscript="/home/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/run_mara_batch_promoters.sh"
runscript="/home/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/run_mara_batch_promoters.sh"
[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

# jmarks="H3K4me3"
# jmarks="H3K4me3"
jmarks="H3K27me3"
# jmarks="H3K4me1"
# jmarks="H3K4me3"

jscale=0
jcenter=0
jbyrow=0

jcenterE="TRUE"
# binarize="TRUE"
# binarize="FALSE"

n=0
maxjobs=1

jmem='96G'
jtime='6:00:00'

jbinskeep=250
jtopics=30
pathprefix="/home/jyeung/hpc"

for jmark in $jmarks; do
    Efname="countmat_PZ_fromHiddenDomains_${jmark}.AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_${jbinskeep}.ntopics_${jtopics}.2020-02-11.txt"
    E="${pathprefix}/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/${jmark}/mara_input/count_mats_peaks_norm/${Efname}"
    Ebase=$(basename $E)
    Ebase=${Ebase%.*}
    [[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

    Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.txt"
    N="${pathprefix}/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/${jmark}/mara_input/sitecount_mats/${Nfname}"
    [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

    Nbase=$(basename $N)
    Nbase=${Nbase%.*}

    outdir="${pathprefix}/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/${jmark}/mara_output/${Ebase}-${Nbase}"
    [[ -d $outdir ]] && echo "$outdir found, continuing for safety" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

    BNAME=$outdir/${jmark}_${dist}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "Running MARA for $jmark"
    echo "bash $runscript $E $outdir $N&"
    bash $runscript $E $outdir $N&
    # echo "bash $runscript $E $outdir $N" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done

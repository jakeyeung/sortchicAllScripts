#!/bin/sh
# Jake Yeung
# 6-run_MARA_on_LDA_output_highK.sh
# Run MARA on LDA output: highK
# 2019-03-13

# # WRAP UP
# while [[ `qstat | grep STDIN | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

runscript="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/run_mara_batch_promoters.sh"
[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

# jmarks="H3K4me3"
# jmarks="H3K4me3"
jmarks="H3K4me1 H3K4me3 H3K27me3"
# jmarks="H3K27me3"
# jmarks="H3K4me3"
# jmarks="H3K4me1"

# keepNbins=0
jscale=0
jcenter=0
jbyrow=0

jcenterE="TRUE"

n=0
maxjobs=1

jmem='32G'
jtime='2:00:00'

# jbinskeep=250
# jtopics=30
# pathprefix="/home/jyeung/hpc"
pathprefix="/hpc/hub_oudenaarden/jyeung/data"
dname="mara_analysis_BM-AllMerged3_Peaks"
marainputdir="${pathprefix}/scChiC/${dname}"

jmem='96G'
jtime='3:00:00'

for jmark in $jmarks; do
    # Efname="BM_${jmark}.BM_AllMerged.VarCorrection.Poisson.GLMPCA.RemoveSmallPeaks.txt"
    Efname="ldaOut.count_mat_from_hiddendomains.${jmark}.filtNAcells_topbins.K-30.keepNbins_0.txt"
    E="${marainputdir}/count_mats_peaks_norm/${Efname}"
    Ebase=$(basename $E)
    Ebase=${Ebase%.*}
    [[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

    # Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.txt"
    # Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.RemoveLargePeaks.${tx}txt"
    Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.${jmark}.txt"
    # RemoveLargePeaks.txt
    N="${marainputdir}/${jmark}/mara_input/sitecount_mats/${Nfname}"
    [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

    Nbase=$(basename $N)
    Nbase=${Nbase%.*}

    outdir="${pathprefix}/scChiC/${dname}/${jmark}/mara_output/${Ebase}-${Nbase}"
    [[ -d $outdir ]] && echo "$outdir found, continuing for safety" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

    BNAME=$outdir/${jmark}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "Running MARA for $jmark"
    echo "bash $runscript $E $outdir $N&"
    # bash $runscript $E $outdir $N&
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; bash $runscript $E $outdir $N" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MARA_${jmark}
    # echo "bash $runscript $E $outdir $N" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done

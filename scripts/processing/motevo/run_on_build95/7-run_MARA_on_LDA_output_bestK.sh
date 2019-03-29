#!/bin/sh
# Jake Yeung
# 6-run_MARA_on_LDA_output_highK.sh
# Run MARA on LDA output: highK
# 2019-03-13

runscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/run_mara_batch_promoters.sh"
[[ ! -e $runscript ]] && echo "$runscript not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmarks="H3K4me3"

jthres="0.99"
jscale=0
jcenter=0
jbyrow=0
jsuffix=""
jsuffixE="Kbest"

jcenterE="TRUE"
binarize="FALSE"  # probably binarize is better for peak? 
# binarize="FALSE"

expersuffix="build95.cells_from_bin_analysis"

jcellmin=0
jcellmax=9999999

dists="0 1000 10000"

n=0
maxjobs=3

jmem='16G'
jtime='1:00:00'

for dist in $dists; do
    for jmark in $jmarks; do
        Efname="hiddenDomains_cellmin_${jcellmin}-cellmax_${jcellmax}-binarize_${binarize}-BM_${jmark}.filt_${jthres}.center_${jcenterE}_${jsuffixE}.txt"
        E="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${expersuffix}/${jmark}/mara_input/count_mats_peaks_norm/$Efname"
        Ebase=$(basename $E)
        Ebase=${Ebase%.*}
        [[ ! -e $E ]] && echo "$E not found, exiting" && exit 1

        Nfname="hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.${dist}.txt"
        N="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${expersuffix}/${jmark}/mara_input/sitecount_mats${jsuffix}/$Nfname"
        [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

        Nbase=$(basename $N)
        Nbase=${Nbase%.*}

        outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${expersuffix}/${jmark}/mara_output/$Ebase-$Nbase-$jsuffix-$jsuffixE"

        [[ -d $outdir ]] && echo "$outdir found, continuing for safety" && continue

        [[ ! -d $outdir ]] && mkdir -p $outdir

        # [[ -d $outdir ]] && echo "$outdir found, exiting to prevent overwrite" && exit 1
        [[ ! -e $N ]] && echo "$N not found, exiting" && exit 1

        BNAME=$outdir/${jmark}_${dist}_qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        # bash $runscript $E $outdir $N&
        echo "bash $runscript $E $outdir $N" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done

#!/bin/sh
# Jake Yeung
# 3-run.computeMatrix.sh
#  
# 2020-02-24

jmem='16G'
jtime='1:00:00'


# jmark="H3K27me3"
# jmark="H3K4me1"
# dist="3000"
dists="3000 5000 10000"
# dists="100000 500000"
jmarks="H3K4me1 H3K4me3 H3K27me3"

jsuffix="FromTopics.2000"
ncores=1
bsizeout=100
bsize="bsize_${bsizeout}"
# bsizeout="10"
basedir="/hpc/hub_oudenaarden/jyeung/data/scChiC"

normby="DefaultNorm"

# dirsuffix=".Downsamp"
dirsuffix=""

refdir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBM${jsuffix}"
[[ ! -d $refdir ]] && echo "$refdir not found, exiting" && exit 1

for jmark in $jmarks; do
    for dist in $dists; do
        echo $jmark

        indir="${basedir}/raw_data/bigwigs_with_R1only/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_40.${bsize}.2020-06-17.offset.r1only${dirsuffix}"
        [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
        # outdir="${basedir}/bigwig_outputs/merged_bams.deeptools_outputs.tss.OneDirPretty.AllGsets.2020-06-17.offset.${jsuffix}.r1only"
        outdir="${basedir}/bigwig_outputs/r1onlyFewerClusters/merged_bams.deeptools_outputs.tss.OneDirPretty.AllGsets.2020-06-17.offset.${jsuffix}.r1only${dirsuffix}"
        [[ ! -d $outdir ]] && mkdir $outdir

        # refs="${refdir}/MouseBM_TSS_FromTopics.Bcell.bsize_2.bed ${refdir}/MouseBM_TSS_FromTopics.Granu.bsize_2.bed ${refdir}/MouseBM_TSS_FromTopics.Eryth.bsize_2.bed ${refdir}/MouseBM_TSS_FromTopics.HSPCs.bsize_2.bed"
        refs="${refdir}/MouseBM_TSS_FromTopics.HSPCs.bsize_2.bed ${refdir}/MouseBM_TSS_FromTopics.Bcell.bsize_2.bed ${refdir}/MouseBM_TSS_FromTopics.Granu.bsize_2.bed ${refdir}/MouseBM_TSS_FromTopics.Eryth.bsize_2.bed"

        # infsinput="${indir}/${jmark}-BM_AllMerged.HSPCs.sorted.${bsizeout}.CPM.bw ${indir}/${jmark}-BM_AllMerged.Bcells.sorted.${bsizeout}.CPM.bw ${indir}/${jmark}-BM_AllMerged.Granulocytes.sorted.${bsizeout}.CPM.bw ${indir}/${jmark}-BM_AllMerged.Erythroblasts.sorted.${bsizeout}.CPM.bw"
        infsinput="${indir}/${jmark}-BM_AllMerged.HSPCs.sorted.cleaned.${bsizeout}.${normby}.bw ${indir}/${jmark}-BM_AllMerged.Bcells.sorted.cleaned.${bsizeout}.${normby}.bw ${indir}/${jmark}-BM_AllMerged.Granulocytes.sorted.cleaned.${bsizeout}.${normby}.bw ${indir}/${jmark}-BM_AllMerged.Erythroblasts.sorted.cleaned.${bsizeout}.${normby}.bw"

        for inf in $infsinput; do
            [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        done

        outf="${outdir}/computeMatrix.${jmark}.${bsize}.${jsuffix}.dist_${dist}.${normby}.tab.gz"
        BNAME=${outf%.*}.AllGsets.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        [[ -e $outf ]] && echo "$outf found, continuing" && continue
        # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; computeMatrix reference-point -S $infsinput -R $refs -bs ${bsizeout} -b ${dist} -a ${dist} -out $outf --skipZeros --missingDataAsZero -p $ncores --sortRegions keep --samplesLabel Bcells Granu Eryth HSPCs"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; computeMatrix reference-point -S $infsinput -R $refs -bs ${bsizeout} -b ${dist} -a ${dist} -out $outf --skipZeros --missingDataAsZero -p $ncores --sortRegions keep --samplesLabel HSPCs Bcells Granu Eryth"
        # --regionsLabel BcellsSpec GranuSpec ErythSpec HSPCsSpec
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"

    done
done


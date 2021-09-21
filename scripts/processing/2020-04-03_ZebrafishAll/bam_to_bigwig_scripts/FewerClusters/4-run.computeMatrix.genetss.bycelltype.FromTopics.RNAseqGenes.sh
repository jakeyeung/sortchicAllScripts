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
# dists="25000 100000"
jmarks="H3K4me1 H3K4me3 H3K27me3"

jsuffix="FromTopics.500"
ncores=4
bsizeout=100
bsize="bsize_${bsizeout}"
# bsizeout="10"
basedir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic"

refdir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/ZebrafishWKM${jsuffix}"
[[ ! -d $refdir ]] && echo "$refdir not found, exiting" && exit 1

for jmark in $jmarks; do
    for dist in $dists; do
        echo $jmark

        # indir="${basedir}/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_40.${bsize}.2020-06-17.offset.r1only"
        indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bigwigs_all/ZebrafishWKM.bigwigs_by_fewerclusters.mapq_40.bsize_100.2020-06-17.offset.r1only"
        [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
        # outdir="${basedir}/bigwig_outputs/merged_bams.deeptools_outputs.tss.OneDirPretty.AllGsets.2020-06-17.offset.${jsuffix}.r1only"
        # outdir="${basedir}/bigwig_outputs/r1onlyFewerClusters/merged_bams.deeptools_outputs.tss.OneDirPretty.AllGsets.2020-06-17.offset.${jsuffix}.r1only"
        outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/deeptoolsout/ZebrafishWKM.deeptools_outputs.2020-06-17.offset.${jsuffix}.r1only"
        [[ ! -d $outdir ]] && mkdir $outdir

        # refs="${refdir}/ZebrafishWKM_TSS_FromTopics.Bcell.bsize_2.bed ${refdir}/ZebrafishWKM_TSS_FromTopics.Granu.bsize_2.bed ${refdir}/ZebrafishWKM_TSS_FromTopics.Eryth.bsize_2.bed ${refdir}/ZebrafishWKM_TSS_FromTopics.HSPCs.bsize_2.bed"
        # ZebrafishWKM_TSS_FromTopics.HSPCs.bsize_2.bed

        refs="${refdir}/ZebrafishWKM_TSS_FromTopics.HSPCs.bsize_2.bed ${refdir}/ZebrafishWKM_TSS_FromTopics.Bcells.bsize_2.bed ${refdir}/ZebrafishWKM_TSS_FromTopics.Granus.bsize_2.bed ${refdir}/ZebrafishWKM_TSS_FromTopics.Eryths.bsize_2.bed"

        for ref in $refs; do
            [[ ! -e $ref ]] && echo "$ref not found, exiting" && exit 1
        done

        # PZ-ChIC-ZF_H3K4me3_2020-04-07.HSPCs.sorted.cleaned.100.cleaned.CPM.bw
        # infsinput="${indir}/PZ-ChIC-ZF_${jmark}_2020-04-07.HSPCs.sorted.cleaned.${bsizeout}.CPM.bw ${indir}/${jmark}-WKM_AllMerged.Bcells.sorted.cleaned.${bsizeout}.CPM.bw ${indir}/${jmark}-WKM_AllMerged.Granulocytes.sorted.cleaned.${bsizeout}.CPM.bw ${indir}/${jmark}-WKM_AllMerged.Erythroblasts.sorted.cleaned.${bsizeout}.CPM.bw"

        # PZ-ChIC-ZF_H3K4me1_2020-04-07.HSPCs.sorted.cleaned.100.cleaned.CPM.bw
        infsinput="${indir}/PZ-ChIC-ZF_${jmark}_2020-04-07.HSPCs.sorted.cleaned.${bsizeout}.cleaned.CPM.bw ${indir}/PZ-ChIC-ZF_${jmark}_2020-04-07.lymph.sorted.cleaned.${bsizeout}.cleaned.CPM.bw ${indir}/PZ-ChIC-ZF_${jmark}_2020-04-07.granu.sorted.cleaned.${bsizeout}.cleaned.CPM.bw ${indir}/PZ-ChIC-ZF_${jmark}_2020-04-07.eryth.sorted.cleaned.${bsizeout}.cleaned.CPM.bw"

        for inf in $infsinput; do
            [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        done

        outf="${outdir}/computeMatrix.${jmark}.${bsize}.${jsuffix}.dist_${dist}.tab.gz"
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


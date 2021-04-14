#!/bin/sh
# Jake Yeung
# 4-run.plotHeatmap.sh
#  
# 2020-02-24

# WRAP UP
while [[ `qstat | grep compMatrix | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='8G'
jtime='1:00:00'

# jmark="H3K4me1"
# jmark="H3K4me1"
# jmark="H3K27me3"
# jdist="3000"
# jdists="3000 5000 10000 25000"
# jmarks="H3K4me1 H3K4me3 H3K27me3"
bsize="bsize_100"
# jsuffix="FromTopics.500"
jsuffix="FromRNAseq"

# for jmark in $jmarks; do
#     for jdist in $jdists; do

        # indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.dist_${jdist}.ctype.sorted.${jmark}.${bsize}.${jsuffix}.AllGsets.2020-06-17"
        # indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.dist_${jdist}.ctype.sorted.${jmark}.${bsize}.${jsuffix}.AllGsets.2020-06-17"
        # indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.OneDir.AllGsets.2020-06-17.${jsuffix}"
        indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.OneDir.AllGsets.2020-06-17.offset.${jsuffix}"
        [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
        outdir=${indir}

        for inf in `ls -d $indir/*.tab.gz`; do
            bname=$(basename $inf)
            bname=${bname%.*}
            jstr=$(echo $bname | cut -d"." -f2,3)
            BNAME=${outdir}/plotHeatmap.${bname}.qsub
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

            outf=${outdir}/plotHeatmap.${bname}.mean.png
            # [[ -e $outf ]] && echo "$outf found, continuing" && continue
            echo $jstr

        # --regionsLabel BcellsSpec GranuSpec ErythSpec HSPCsSpec

            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; plotHeatmap -m $inf -out $outf --colorMap YlGnBu --heatmapHeight 15 --sortRegions descend --sortUsingSamples 1 --averageTypeSummaryPlot mean --regionsLabel BcellsSpec GranuSpec ErythSpec HSPCsSpec --samplesLabel Bcells Granu Eryth HSPCs --legendLocation best"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
        done

#     done
# done



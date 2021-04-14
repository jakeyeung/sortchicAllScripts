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

jmem='4G'
jtime='0:30:00'

markprefix="refmark_H3K27me3"

# jmark="H3K4me1"
# jmark="H3K4me1"
# jmark="H3K27me3"
# jdist="3000"
# jdists="3000 5000 10000 25000"
# jmarks="H3K4me1 H3K4me3 H3K27me3"
bsize="bsize_100"
jsuffix="FromTopics.${markprefix}.2000"
# dirsuffix=".Downsamp"
dirsuffix=""
# dirsuffix2=".NormHeteroCuts"
# dirsuffix2=".NormHeteroCutsNcells"
# dirsuffix2=".NormTotalCuts"
# dirsuffix2=".NormTotalCutsOnly"
# dirsuffix2=".NormHeteroCutsOnly"
# dirsuffix2=".NormHeteroCutsPerTotalOnly"
# dirsuffix2=".NormTotalCutsPerCell"
# dirsuffix2=".NormHeteroCutsPerCell"

# dirsuffix2=".MillionBinCuts"
# dirsuffix2=".MillionTSSCuts"
# dirsuffix2=".TSSCuts"
dirsuffix3=".DefaultNorm"

dirsuffix2s=".MillionTSSCuts .MillionBinCuts"

for dirsuffix2 in $dirsuffix2s; do
        indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/r1onlyFewerClusters/merged_bams.deeptools_outputs.tss.OneDirPretty.AllGsets.2020-06-17.offset.${jsuffix}.r1only${dirsuffix}${dirsuffix2}${dirsuffix3}.${markprefix}"
        # "merged_bams.deeptools_outputs.tss.OneDirPretty.AllGsets.2020-06-17.offset.FromTopics.refmark_H3K27me3.2000.r1only.MillionBinCuts.DefaultNorm.refmark_H3K27me3"

        [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
        outdir=${indir}

        for inf in `ls -d $indir/*.tab.gz`; do
            bname=$(basename $inf)
            bname=${bname%.*}
            jstr=$(echo $bname | cut -d"." -f2,3)
            BNAME=${outdir}/plotHeatmap.${bname}.qsub
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

            outf=${outdir}/plotHeatmap.${bname}.mean${dirsuffix}.${markprefix}.png
            # [[ -e $outf ]] && echo "$outf found, continuing" && continue
            echo $jstr

            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; plotHeatmap -m $inf -out $outf --colorMap YlGnBu --heatmapHeight 15 --sortRegions descend --sortUsingSamples 1 --averageTypeSummaryPlot mean --regionsLabel HSPCsSpec BcellsSpec GranuSpec ErythSpec --samplesLabel HSPCs Bcells Granu Eryth --legendLocation best --yAxisLabel CountsPer${dirsuffix2} --plotTitle ${jstr}" 
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
        done


done
#     done
# done



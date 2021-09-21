#!/bin/sh
# Jake Yeung
# 2-bam_to_bigwig.sh
#  
# 2020-01-10

# # WRAP UP
# while [[ `qstat | grep split | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='16G'
jtime='6:00:00'

mapq="40"

suffix2=".ScaleByTssOverBinCuts"
# suffix=".Downsamp"
suffix=""
# bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.sh"
bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.WithScaleFactor.CPM.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_${mapq}.2020-06-17"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits${suffix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

# indirsf="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts.countsByCluster"
indirsf="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bins_tss_totalcuts.countsByCluster"

# bsizes="100"
bsize="100"
# for bsize in $bsizes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/bigwigs_with_R1only/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_${mapq}.bsize_${bsize}.2020-06-17.offset.r1only${suffix}${suffix2}.Verbose"

    # [[ ! -d $outdir ]] && mkdir $outdir
    [[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1
    for inbam in `ls -d $indir/*.bam`; do
        echo $inbam
        bname=$(basename $inbam)
        bname=${bname%.*}  # strip extension 
        outbw="${outdir}/${bname}.${bsize}.CPM${suffix2}.bw"

        [[ -e $outbw ]] && echo "$outbw found, continuing" && continue

        BNAME=$outdir/$bname.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        # get SF
        clst=$(echo $bname | cut -d"." -f2)
        jmark=$(echo $bname | cut -d"-" -f1)
        if [ $jmark == "H3K9me3" ]
        then
            echo "Skipping H3K9me3 file: $inbam"
        	continue
        fi  

        # echo "$clst $jmark"
        # infsf="${indirsf}/hetero_and_totalcuts.${clst}.${jmark}.NoCname.txt"
        infsf="${indirsf}/MouseBMFromTopics.2000.bins_tss_hetero_and_totalcuts.${clst}.${jmark}.NoCname.txt"
        [[ ! -e $infsf ]] && echo "$infsf not found, exiting" && exit 1
        bincuts=`awk '{print $5}' $infsf`  # cluster.new     ncuts.total     ncuts.hetero    tss.cuts        bincuts     ncells  mark
        tsscuts=`awk '{print $4}' $infsf`  # cluster.new     ncuts.total     ncuts.hetero    tss.cuts        bincuts    ncells  mark
        ncells=`awk '{print $6}' $infsf`

        sf=`awk -v bincuts=$bincuts -v tsscuts=$tsscuts -v ncells=$ncells 'BEGIN { print (tsscuts / bincuts)  }'`
        # fullsf=`awk -v ncells=$ncells -v bincuts=$bincuts -v tsscuts=$tsscuts 'BEGIN { print ( 1 /  tsscuts ) }'`
        echo "Reading scale factors from $infsf"
        echo "Calculating $tsscuts/$bincuts $sf"

        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N "bam2bw_${bname}"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl $sf"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
    done
# done

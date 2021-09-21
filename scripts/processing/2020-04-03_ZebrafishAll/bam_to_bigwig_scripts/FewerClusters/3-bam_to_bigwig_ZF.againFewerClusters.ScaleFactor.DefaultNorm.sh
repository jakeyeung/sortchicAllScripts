# !/bin/sh
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

normby="DefaultNorm"
# suffix2=".ScaleByNcells"  # 4th column in infsf
# suffix2=".NormHeteroCuts"  # 3rd column in infsf
# suffix2=".NormHeteroCutsNcells"  # 3rd and 4th column in infsf
# suffix2=".NormTotalCuts"  # 2nd column infsf
# suffix2=".NormTotalCutsOnly"  # 2nd column infsf
# suffix2=".NormHeteroCutsOnly"  # 2nd column infsf
# suffix2=".NormHeteroCutsPerTotalOnly"  # 2nd column infsf

# suffix2=".NormTotalCutsPerCell"  # 2nd column * 4th column infsf
# suffix2=".NormHeteroCutsPerCell"  # 2nd column * 4th column infsf
# suffix2=".MillionBinCuts"
suffix2=".MillionTSSCuts"

suffix=""

bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.WithScaleFactor.${normby}.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

jprefix="imputevarfilt.lessstringent"
suffix=""
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17/DedupR1onlyNoAltHits${suffix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/zebrafish/WKM/correlated_bins.chromofixed.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

# indirsf="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts.countsByCluster"
indirsf="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bins_tss_totalcuts.countsByCluster/zebrafish"
[[ ! -d $indirsf ]] && echo "$indirsf not found, exiting" && exit 1

# bsizes="100"
bsize="100"
# for bsize in $bsizes; do
    # outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/bigwigs_with_R1only/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_${mapq}.bsize_${bsize}.2020-06-17.offset.r1only${suffix}${suffix2}.${normby}"
    outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bigwigs_all/ZebrafishWKM.bigwigs_by_fewerclusters.mapq_${mapq}.bsize_${bsize}.2020-06-17.offset.r1only${suffix}${suffix2}"
    [[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1
    # [[ ! -d $outdir ]] && mkdir $outdir

    for inbam in `ls -d $indir/*.bam`; do
        echo $inbam
        bname=$(basename $inbam)
        bname=${bname%.*}  # strip extension 
        outbw="${outdir}/${bname}.${bsize}.${normby}${suffix2}.bw"

        [[ -e $outbw ]] && echo "$outbw found, continuing" && continue

        BNAME=$outdir/$bname.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        # get SF
        clst=$(echo $bname | cut -d"." -f2)
        jmark=$(echo $bname | cut -d"_" -f2)
        if [ $jmark == "H3K9me3" ]
        then
            echo "Skipping H3K9me3 file: $inbam"
        	continue
        fi  

        echo "Cluster: $clst"
        echo "Mark: $jmark"
        # echo "Debug exit"
        # exit 0

        # echo "$clst $jmark"
        # infsf="${indirsf}/MouseBMFromTopics.2000.bins_tss_hetero_and_totalcuts.${clst}.${jmark}.NoCname.txt"
        infsf="${indirsf}/ZebrafishWKMFromTopics.2000.bins_tss_hetero_and_totalcuts.${clst}.${jmark}.NoCname.txt"
        # infsf="${indirsf}/hetero_and_totalcuts.${clst}.${jmark}.NoCname.txt"
        [[ ! -e $infsf ]] && echo "$infsf not found, exiting" && exit 1
        bincuts=`awk '{print $5}' $infsf`  # cluster.new     ncuts.total     ncuts.hetero    tss.cuts        bincuts     ncells  mark
        tsscuts=`awk '{print $4}' $infsf`  # cluster.new     ncuts.total     ncuts.hetero    tss.cuts        bincuts    ncells  mark
        ncells=`awk '{print $6}' $infsf`
        # echo "$infsf $sf"

        # fullsf=`awk -v nf=$nf -v ncells=$ncells -v bincuts=$bincuts 'BEGIN { print ( 1 / bincuts ) }'`
        # fullsf=$(awk -v ncells=$ncells -v tsscuts=$tsscuts 'BEGIN { print ( ncells / tsscuts ) }')

        echo "Reading scale factors from: $infsf"

        echo "ncells: $ncells bincuts: $bincuts tsscuts: $tsscuts"
        fullsf=`awk -v ncells=$ncells -v bincuts=$bincuts -v tsscuts=$tsscuts 'BEGIN { print ( 1000000 /  tsscuts ) }'`
        # fullsf=`awk -v ncells=$ncells -v bincuts=$bincuts -v tsscuts=$tsscuts 'BEGIN { print ( 1000000 /  bincuts ) }'`
        echo "Scale factor: $fullsf"
        # echo "Debug exit"
        # exit 0
        
        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N "bam2bw_${bname}"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl $fullsf"
        # $cmd
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
        # exit 0
    done
# done

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
suffix2=".NormHeteroCutsPerCell"  # 2nd column * 4th column infsf
# suffix2="NormTSSCuts"

# suffix=".Downsamp"
suffix=""
# bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.sh"
bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.WithScaleFactor.${normby}.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_${mapq}.2020-06-17"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits${suffix}"
# indirDefaultNorm="$indir/default_norm_factors"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

indirsf="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts.countsByCluster"

# bsizes="100"
bsize="100"
# for bsize in $bsizes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/bigwigs_with_R1only/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_${mapq}.bsize_${bsize}.2020-06-17.offset.r1only${suffix}${suffix2}.${normby}"
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

        # get default normfac

        # infnf="$indirDefaultNorm/${bname}.DefaultNormFactor.txt"
        # [[ ! -e $infnf ]] && echo "$infnf not found, exiting" && exit 1
        # nf=`awk '{print $1}' $infnf`
        # # echo $nf

        # get SF
        clst=$(echo $bname | cut -d"." -f2)
        jmark=$(echo $bname | cut -d"-" -f1)
        if [ $jmark == "H3K9me3" ]
        then
            echo "Skipping H3K9me3 file: $inbam"
        	continue
        fi  

        # echo "$clst $jmark"
        infsf="${indirsf}/hetero_and_totalcuts.${clst}.${jmark}.NoCname.txt"
        [[ ! -e $infsf ]] && echo "$infsf not found, exiting" && exit 1
        totalcuts=`awk '{print $2}' $infsf`  # colnames: ctype, totalcuts, heterocuts, ncells, mark
        # sf=`awk '{print $3}' $infsf`  # colnames: ctype, totalcuts, heterocuts, ncells, mark
        heterocuts=`awk '{print $3}' $infsf`  # colnames: ctype, totalcuts, heterocuts, ncells, mark
        ncells=`awk '{print $4}' $infsf`  # colnames: ctype, totalcuts, heterocuts, ncells, mark
        # echo "$infsf $sf"

        # full scale factor
        # echo "$nf / ($sf / $ncells)"
        # fullsf=`awk -v nf=$nf -v sf=$sf -v ncells=$ncells 'BEGIN { print ( nf / ( sf / ncells ) ) }'`
        # echo $fullsf

        # fullsf=`awk -v nf=$nf -v totalcuts=$totalcuts 'BEGIN { print ( 1 / totalcuts ) }'`
        # fullsf=`awk -v nf=$nf -v heterocuts=$heterocuts 'BEGIN { print ( 1 / heterocuts ) }'`
        # fullsf=`awk -v nf=$nf -v heterocuts=$heterocuts -v totalcuts=$totalcuts 'BEGIN { print ( totalcuts / heterocuts ) }'`
        # fullsf=`awk -v nf=$nf -v ncells=$ncells -v totalcuts=$totalcuts 'BEGIN { print ( ncells / totalcuts ) }'`
        # fullsf=`awk -v nf=$nf -v ncells=$ncells -v heterocuts=$heterocuts 'BEGIN { print ( ncells / heterocuts ) }'`
        fullsf=`awk -v nf=$nf -v ncells=$ncells -v bincuts=$bincut 'BEGIN { print ( bincuts ) }'`
        echo $fullsf

        # echo "Debug. exiting"
        # exit 0

        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N "bam2bw_${bname}"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl $fullsf"
        # $cmd
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
    done
# done

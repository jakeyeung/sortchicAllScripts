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
dists="3000 5000 10000 25000"
jmarks="H3K4me1 H3K4me3 H3K27me3"

jsuffix="FromTopics.500"
ncores=4
bsizeout=10
bsize="bsize_${bsizeout}"
# bsizeout="10"
basedir="/hpc/hub_oudenaarden/jyeung/data/scChiC"

refdir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBM${jsuffix}"
[[ ! -d $refdir ]] && echo "$refdir not found, exiting" && exit 1

for jmark in $jmarks; do
    for dist in $dists; do
        echo $jmark

        indir="${basedir}/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_40.${bsize}.2020-06-17"
        # outdir="${basedir}/bigwig_outputs/merged_bams.deeptools_outputs.tss.dist_${dist}.ctype.sorted.${jmark}.${bsize}.${bsizeout}.${jsuffix}.2020-06-17"
        outdir="${basedir}/bigwig_outputs/merged_bams.deeptools_outputs.tss.dist_${dist}.ctype.sorted.${jmark}.${bsize}.${jsuffix}.2020-06-17"
        [[ ! -d $outdir ]] && mkdir $outdir

        for ref in `ls -d $refdir/*.bed`; do
            bname=$(basename $ref)
            bname=${bname%.*}
            ctype=$(echo $bname | cut -d"." -f2)
            infsinput=${indir}/${jmark}*.bw

            if [ $ctype == "LowExprs" ]
            then
                echo "Skipping $ctype"
                continue
            fi

            echo $ctype
            # continue

            outf="${outdir}/computeMatrix.${jmark}.${ctype}.${jsuffix}.sorted.tab.gz"
            BNAME=${outf%.*}.${ctype}.qsub
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
            [[ -e $outf ]] && echo "$outf found, continuing" && continue
            cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; computeMatrix reference-point -S $infsinput -R $ref -bs ${bsizeout} -b ${dist} -a ${dist} -out $outf --skipZeros --missingDataAsZero -p $ncores --sortRegions keep"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
        done
    done
done


#!/bin/sh
# Jake Yeung
# 3-run.computeMatrix.sh
#  
# 2020-02-24

jmem='16G'
jtime='1:00:00'


jmark="H3K27me3"
jsuffix="genetss"

ncores=4

# jsuffix="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.bsize_500"
bsize="bsize_100"
dist="5000"
basedir="/hpc/hub_oudenaarden/jyeung/data/scChiC"
indir="${basedir}/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.${bsize}"
outdir="${basedir}/bigwig_outputs/merged_bams.deeptools_outputs.tss.dist_${dist}.ctype.sorted.${jmark}.${bsize}"
[[ ! -d $outdir ]] && mkdir $outdir

# refdir="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered_sorted"
refdir="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered.ctypes"

for ref in `ls -d $refdir/*.bed`; do
    bname=$(basename $ref)
    bname=${bname%.*}
    ctype=$(echo $bname | cut -d"." -f5 | cut -d"_" -f2)
    infsinput=${indir}/${jmark}*.bw

    echo $ctype
    # continue

    outf="${outdir}/computeMatrix.${jmark}.${ctype}.${jsuffix}.sorted.tab.gz"
    BNAME=${outf%.*}.${ctype}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; computeMatrix reference-point -S $infsinput -R $ref -bs 100 -b ${dist} -a ${dist} -out $outf --skipZeros --missingDataAsZero -p $ncores --sortRegions keep" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N compMatrix_${jmark}_${jsuffix}
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; computeMatrix reference-point -S $infsinput -R $ref -bs 100 -b ${dist} -a ${dist} -out $outf --skipZeros --missingDataAsZero -p $ncores --sortRegions keep"
done

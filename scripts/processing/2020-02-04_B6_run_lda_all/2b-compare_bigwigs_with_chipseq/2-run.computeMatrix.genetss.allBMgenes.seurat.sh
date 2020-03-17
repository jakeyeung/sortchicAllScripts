#!/bin/sh
# Jake Yeung
# 3-run.computeMatrix.sh
#  
# 2020-02-24

jmem='16G'
jtime='1:00:00'


# jmark="H3K27me3"
# jmark="H3K4me1"
jmark="H3K4me3"

ncores=4

bs=100
bsize="bsize_100"
dist="5000"
basedir="/hpc/hub_oudenaarden/jyeung/data/scChiC"
mapq="40"
# indir="${basedir}/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.${bsize}"
indir="${basedir}/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.MAPQ_${mapq}.${bsize}"
outdir="${basedir}/bigwig_outputs/merged_bams.deeptools_outputs.tss.MAPQ_${mapq}.dist_${dist}.allctypes_from_seurat.${jmark}.${bsize}"
[[ ! -d $outdir ]] && mkdir $outdir

refdir="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered.seurat"

for ref in $(ls -d $refdir/*.bed); do
    echo $ref
    jbase=$(basename $ref)
    jbase=${jbase%.*}

    bname=$(basename $ref)
    bname=${bname%.*}
    ctype=$(echo $bname | cut -d"." -f5 | cut -d"_" -f2)
    infsinput=${indir}/${jmark}*.bw

    outf="${outdir}/computeMatrix.MAPQ_${mapq}.${jmark}.${jbase}.tab.gz"
    BNAME=${outf%.*}.allctypes.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; computeMatrix reference-point -S $infsinput -R $ref -bs ${bs} -b ${dist} -a ${dist} -out $outf --skipZeros --missingDataAsZero -p $ncores --sortRegions keep" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N compMatrix_${jmark}_${jbase}


done

# ref="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/giladi_filtered.seurat/gene_tss_winsize.2.BM_giladi.all_genes_from_seurat.pval_0.01.logfc_0.high.bed"
# [[ ! -e $ref ]] && echo "$ref not found, exiting" && exit 1

#!/bin/sh
# Jake Yeung
# 0-remove_eryth_from_mats.sh
#  
# 2020-08-25

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/hiddendomains_analysis/remove_cluster_from_mat.R"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains.NewCountFilters.rds_format"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir="${indir}/remove_eryth"
[[ ! -d $outdir ]] && mkdir $outdir

annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables"
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    inbase="merged.${jmark}.minlength_1000.cutoff_analysis.merged.withchr.annotated.NewCountFilters.rds"
    inf=${indir}/${inbase}
    bname=${inbase%.*}
    outfile=${outdir}/${bname}.remove_eryth.rds
    [[ -e $outfile ]] && echo "$outfile found, continuing" && continue
    annotfile="${annotdir}/BM_AllMerged.${jmark}.cell_cluster_table.txt"
    [[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -annotfile $annotfile -outfile $outfile -clstremove Eryth" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -annotfile $annotfile -outfile $outfile -clstremove Eryth
    # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -annotfile $annotfile -outfile $outfile -clstremove Eryth"
    # echo $cmd
done

#!/bin/sh
# Jake Yeung
# 2-run.merge_count_tables_into_one.sh
#  
# 2021-02-15

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/dynamic_bins_analysis/make_bins/merge_count_tables_into_one.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/count_tables_all_four_marks_dynamic_bins"
outdir="$indir/merged_by_mark"
indirbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.dynamic_bins_TSS_TES_regions"

[[ ! -d $outdir ]] && mkdir $outdir
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# merge TSS bins 

cd $indirbed

metadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"

for mark in $marks; do
    metaname="metadata_batch_corrected.arranged_by_lineage.${mark}.2021-01-02.txt"
    metafile=${metadir}/${metaname}
    [[ ! -e $metafile ]] && echo "$metafile not found, exiting" && exit 1
    # bedfnames="dynamic_bins.50kb.${mark}.txt celltype_specific_genes.TSS_TES.txt celltype_specific_genes.TSS_10kb.txt"
    bedfnames="celltype_specific_genes.TSS_10kb.txt"
    for bedname in $bedfnames; do
        echo $bedname
        infs=${indir}/*${mark}*${bedname}
        echo $infs
        outf="${outdir}/count_tables_merged.${mark}.${bedname}.rds"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $infs -outfile $outf -metafile $metafile"
        echo $cmd
        BNAME=${outdir}/sbatch_${mark}_${bedname}
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${mark}_${bedname} --wrap "$cmd"
    done
done


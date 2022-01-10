#!/bin/sh
# Jake Yeung
# copy_genomewide_counts_for_scchix.sh
#  
# 2022-01-03

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/from_genomewide"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/scchix_inputs_H34me1_H3K9me3_genomewide"

jmarks="H3K4me1 H3K9me3 H3K4me1-H3K9me3"

for jmark in $jmarks; do
    inf="${indir}/count_name.${jmark}.k4_k9_50kb_genomewide.2021-01-30.rds"
    outf="${outdir}/countmat_var_filt.${jmark}.rds"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1
    cp $inf $outf
done


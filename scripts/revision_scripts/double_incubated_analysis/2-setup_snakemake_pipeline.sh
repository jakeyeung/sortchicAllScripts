#!/bin/sh
# Jake Yeung
# 2-setup_snakemake_pipeline.sh
#  
# 2022-01-03

template="/hpc/hub_oudenaarden/jyeung/data/dblchic/es_npc/snakemake/snakemake_pipeline_common_rows_day24only_binfilt_cellfilt/"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/snakemake_H3K4me1_H3K9me3"

# rsync -avrL --copy-links --include "*/"  --include="*.sh" --exclude="*" $template $outdir
# rsync -avrL --copy-links --include "*/"  --include="*.yaml" --exclude="*" $template $outdir
# rsync -avrL --copy-links --include "*/"  --include="Snakefile" --exclude="*" $template $outdir
# cp $template/Snakefile $outdir/Snakefile
cp $template/cluster.json $outdir/cluster.json


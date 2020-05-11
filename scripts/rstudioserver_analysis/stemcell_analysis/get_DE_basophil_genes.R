# Jake Yeung
# Date of Creation: 2020-05-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/get_DE_basophil_genes.R
# Get Basophil for pEter


jsuffix <- ".celltypes_filt"

pval.min <- 0.01
fc.min <- 0

inf.de <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat", jsuffix, ".rds")
assertthat::assert_that(file.exists(inf.de))
dat.de <- readRDS(inf.de)

baso.genes <- subset(dat.de, cluster == "Prss34" & p_val_adj <= pval.min & avg_logFC > fc.min)

# write to output

outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/genesets/basophil_specific_genes_seurat_table.txt"
fwrite(baso.genes, file = outf, sep = "\t")

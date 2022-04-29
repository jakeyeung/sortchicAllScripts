# Jake Yeung
# Date of Creation: 2022-04-15
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/21-plot_expression_on_umaps_clean_eryths.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)


# Load umap  --------------------------------------------------------------

# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_keep_eryths/metadata_celltyping_k27me3.dynamicbins.2022-04-13.txt"
inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.k27me3.2022-04-15.txt"
dat.meta <- fread(inf.meta)


# Get eryth cluster -------------------------------------------------------

dat.meta.eryths <- subset(dat.meta, ctype.from.LL == "Eryths")

# ggplot(dat.meta.eryths, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
#   geom_point() + 
#   scale_color_viridis_c(direction = -1) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 


# Get TSS imputed ---------------------------------------------------------

# inf.lda <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_keep_eryths/lda_outputs.count_mat_merged_with_old_TSS.k27me3.2022-04-12/ldaOut.count_mat_merged_with_old_TSS.k27me3.2022-04-12.Robj"
inf.lda <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_TSS.k27me3.2022-04-12/ldaOut.count_mat_merged_with_old_TSS.k27me3.2022-04-12.Robj"
load(inf.lda, v=T)

tm <- AddTopicToTmResult(posterior(out.lda))

dat.impute <- t(log2(tm$topics %*% tm$terms))

# Plot marker genes -------------------------------------------------------

# get marker genes for eryths and plot

jgene <- "Sox6"
jgene <- "Hbb"
jgene <- "Gata1"
# jgene <- "Cd34"

markergenes <- c("Hbb-bt", "Sox6", "Tal1", "Gata1", "Ebf1", "Cd79a", "Cd79b", "Hoxa9", "Meis1", "Runx2", "Kit", "Hlf", "Erg$", "Cd34", "Stat4", "Tcf7", "Cebpe", "Ly6c2", "Ly6g", "S100a8", "S100a2", "Chil3")

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_no_frip_filter_on_eryths_clean2/BM_k27me3"
dir.create(outdir, recursive = TRUE)

# outdirpdf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_clean_eryths"
# dir.create(outdirpdf)
outpdf <- file.path(outdir, paste0("check_marker_genes_on_UMAP_clean_eryths.", Sys.Date(), ".pdf"))
pdf(outpdf, useDingbats = FALSE)



ggplot(dat.meta, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.eryths, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.eryths, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




for (jgene in markergenes){
  (jcoords <- rownames(dat.impute)[which(grepl(jgene, rownames(dat.impute)))])
  jcoord <- jcoords[[1]]
  jrow <- dat.impute[jcoord, ]
  
  dat.exprs <- data.frame(cell = names(jrow), signal = jrow, stringsAsFactors = FALSE) %>%
    left_join(dat.meta, ., by = "cell")
  dat.exprs.sub <- data.frame(cell = names(jrow), signal = jrow, stringsAsFactors = FALSE) %>%
    left_join(dat.meta.eryths, ., by = "cell")
  
  m <- ggplot(dat.exprs, aes(x = umap1, y = umap2, color = signal)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    ggtitle(jcoord) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m.sub <- ggplot(dat.exprs.sub, aes(x = umap1, y = umap2, color = signal)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c() + 
    ggtitle(jcoord) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.sub)
  
  # m.sub.dens <- ggplot(dat.exprs.sub, aes(x = signal)) + 
  #   geom_density() + 
  #   theme_bw() + 
  #   scale_color_viridis_c() + 
  #   ggtitle(jcoord) + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # print(m.sub.dens)
  
}


# Remove bad eryth cluster  -----------------------------------------------
# 
# m.check <- subset(dat.meta.eryths, batch == "New") %>%
#   ggplot(., aes(x = cell.var.within.sum.norm)) + 
#   geom_density() 
# print(m.check)
# 
# # dat.meta.eryths.label <- dat.meta.eryths %>%
# #   mutate(highlight = batch == "New" & ctype == "Eryths" & cell.var.within.sum.norm < 1.5)
# 
# dat.meta.eryths.label <- dat.meta.eryths %>%
#   mutate(is.bad = batch == "New" & cell.var.within.sum.norm < 1.25 & umap1 > 5)
# 
# ggplot(dat.meta.eryths.label, aes(x = umap1, y = umap2, color = is.bad)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




dev.off()


# Create bad cells  -------------------------------------------------------


bad.eryths <- subset(dat.meta, ctype.from.LL == "Eryths" & umap1 < 20)$cell
bad.basos <- subset(dat.meta, ctype.from.LL == "Basophils" & umap2 < 12)$cell
bad.bcells <- subset(dat.meta, ctype.from.LL == "Bcells" & umap1 > -18)$cell
# bad.dcs <- subset(dat.meta, ctype.from.LL == "DCs" & umap2 < 12)$cell
bad.granus <- subset(dat.meta, ctype.from.LL == "Granulocytes" & umap2 < 12)$cell
bad.hscs <- subset(dat.meta, ctype.from.LL == "HSCs" & umap2 > 0)$cell
bad.meps <- subset(dat.meta, ctype.from.LL == "MEP" & umap1 < 2)$cell
bad.nks <- subset(dat.meta, ctype.from.LL == "NKs" & umap2 > -25)$cell

# check each
bad.cells.lst <- list("Eryths" = bad.eryths, 
                      "Basophils" = bad.basos, 
                      "Bcells" = bad.bcells, 
                      # "DCs" = bad.dcs, 
                      "Granulocytes" = bad.granus, 
                      "HSCs" = bad.hscs, 
                      "MEP" = bad.meps, 
                      "NKs" = bad.nks)

bad.names <- names(bad.cells.lst); names(bad.names) <- bad.names

m.check <- lapply(bad.names, function(jname){
  bad.cells <- bad.cells.lst[[jname]]
  m <- ggplot(dat.meta %>% mutate(is.bad = cell %in% bad.cells), aes(x = umap1, y = umap2, color = is.bad)) + 
    geom_point() + 
    ggtitle(jname) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
})

bad.cells <- unlist(bad.cells.lst)

cells.keep <- dat.meta$cell[!dat.meta$cell %in% bad.cells]

dat.meta.new <- subset(dat.meta, batch == "New")
cells.keep.new <- dat.meta.new$cell[!dat.meta.new$cell %in% bad.cells]

print(dim(dat.meta))
print(length(cells.keep))
print(length(bad.cells))

ggplot(dat.meta %>% mutate(is.bad = cell %in% bad.cells), aes(x = umap1, y = umap2, color = is.bad)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())





# Load and save mats filter out bad cells ---------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_no_frip_filter_on_eryths_clean/BM_k27me3"

saveRDS(bad.cells.lst, file.path(outdir, paste0("bad_cells_manually_removed.k27me3.", Sys.Date(), ".rds")))
fwrite(dat.meta, file.path(outdir, paste0("metadata.k27me3.", Sys.Date(), ".txt")))


inf.allcells.allbins <- file.path(indir, "count_mat_merged_with_old_allbins.k27me3.2022-04-12.rds")
inf.newcells.allbins <- file.path(indir, "count_mat_new_only_allbins.k27me3.2022-04-12.rds")

inf.allcells.dynbins <- file.path(indir, "count_mat_merged_with_old_dynbins.k27me3.2022-04-12.rds")
inf.newcells.dynbins <- file.path(indir, "count_mat_new_only_dynbins.k27me3.2022-04-12.rds")

inf.allcells.tss <- file.path(indir, "count_mat_merged_with_old_TSS.k27me3.2022-04-12.rds")
inf.newcells.tss <- file.path(indir, "count_mat_new_only_TSS.k27me3.2022-04-12.rds")

jdateout <- Sys.Date()
outf.allcells.allbins <- file.path(outdir, paste0("count_mat_merged_with_old_allbins.k27me3.", jdateout, ".rds"))
outf.newcells.allbins <- file.path(outdir, paste0("count_mat_new_only_allbins.k27me3.", jdateout, ".rds"))
outf.allcells.dynbins <- file.path(outdir, paste0("count_mat_merged_with_old_dynbins.k27me3.", jdateout, ".rds"))
outf.newcells.dynbins <- file.path(outdir, paste0("count_mat_new_only_dynbins.k27me3.", jdateout, ".rds"))
outf.allcells.tss <- file.path(outdir, paste0("count_mat_merged_with_old_TSS.k27me3.", jdateout, ".rds"))
outf.newcells.tss <- file.path(outdir, paste0("count_mat_new_only_TSS.k27me3.", jdateout, ".rds"))


mat.allcells.allbins <- readRDS(inf.allcells.allbins)
mat.newcells.allbins <- readRDS(inf.newcells.allbins)

mat.allcells.dynbins <- readRDS(inf.allcells.dynbins)
mat.newcells.dynbins <- readRDS(inf.newcells.dynbins)

mat.allcells.tss <- readRDS(inf.allcells.tss)
mat.newcells.tss <- readRDS(inf.newcells.tss)

mat.allcells.allbins.filt <- mat.allcells.allbins[, cells.keep]
mat.newcells.allbins.filt <- mat.newcells.allbins[, cells.keep.new]

mat.allcells.dynbins.filt <- mat.allcells.dynbins[, cells.keep]
mat.newcells.dynbins.filt <- mat.newcells.dynbins[, cells.keep.new]

mat.allcells.tss.filt <- mat.allcells.tss[, cells.keep]
mat.newcells.tss.filt <- mat.newcells.tss[, cells.keep.new]

print(mat.allcells.allbins.filt)

saveRDS(mat.allcells.allbins.filt, file = outf.allcells.allbins)
saveRDS(mat.newcells.allbins.filt, file = outf.newcells.allbins)

saveRDS(mat.allcells.dynbins.filt, file = outf.allcells.dynbins)
saveRDS(mat.newcells.dynbins.filt, file = outf.newcells.dynbins)

saveRDS(mat.allcells.tss.filt, file = outf.allcells.tss)
saveRDS(mat.newcells.tss.filt, file = outf.newcells.tss)





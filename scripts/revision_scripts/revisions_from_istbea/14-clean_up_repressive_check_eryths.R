# Jake Yeung
# Date of Creation: 2022-04-05
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/14-clean_up_repressive_check_eryths.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k27me3", "k9me3"); names(jmarks) <- jmarks

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_repressive_cleaned_check_eryths_keep_old"
dir.create(outdir)
outpdf <- file.path(outdir, paste("repressive_marks_cleanup_check_eryths_k27me3.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# Load repressive umaps with some bad cells -------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots"

indir.var <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k27me3_var_filt"

infs.meta <- lapply(jmarks, function(jmark){
  print(jmark)
  file.path(indir, paste0("metadata_with_colors.", jmark, ".2022-02-13.txt"))
})

dat.meta.colors <- fread("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap/metadata.k4me1.2022-03-01.txt") %>%
  dplyr::select(ctype.from.LL, colcode) 
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

infs.meta.var <- lapply(jmarks, function(jmark){
  file.path(indir.var, paste0("metadata_with_varfilt2.", jmark, ".txt"))
})

dat.meta1.lst <- lapply(infs.meta, function(jinf) fread(jinf))

dat.meta.var.lst <- lapply(infs.meta.var, function(jinf) fread(jinf) %>% dplyr::select(cell, cell.var.within.sum.norm))

dat.meta.lst <- lapply(jmarks, function(jmark){
  left_join(dat.meta1.lst[[jmark]], dat.meta.var.lst[[jmark]])
})
  
lapply(dat.meta.lst, dim)

# Remove bad cells --------------------------------------------------------

ggplot(dat.meta.lst$k27me3, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k27me3") + 
  scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                        guide = "legend") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#3adf31", "#f3f119", "#e1383e", "#5c05f5", "#aab15e", "#b2943b")
ggplot(dat.meta.lst$k27me3, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k27me3") + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.lst$k27me3, aes(x = umap1, y = umap2, color = ctype)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~ctype) + 
  ggtitle("k27me3") + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.lst$k27me3 %>% mutate(highlight = ctype.from.LL == "NKs"), 
       aes(x = umap1, y = umap2, color = highlight)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~batch) + 
  geom_vline(xintercept = c(13, 23)) + 
  ggtitle("k27me3") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.lst$k27me3, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k27me3") + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#99a339")
ggplot(dat.meta.lst$k27me3, aes(x = umap1, y = umap2, color = as.character(louvain))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k27me3") + 
  scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Show umap cleaned up  ---------------------------------------------------

# bad.louvs.k27 <- c("16", "12")
bad.louvs.k27 <- c("16", "12")
# bad.louvs.k27 <- c()

jsub.k27 <- dat.meta.lst$k27me3 %>%
  mutate(is.bad = louvain %in% bad.louvs.k27 & batch == "New")

# ggplot(dat.meta.lst$k27me3 %>%
#          mutate(is.bad = louvain %in% bad.louvs.k27 & batch == "New"), 
ggplot(jsub.k27, 
       aes(x = umap1, y = umap2, color = is.bad)) + 
  geom_point() + 
  facet_wrap(~batch) + 
  theme_bw() + 
  geom_vline(xintercept = 15) +
  ggtitle("k27me3") + 
  scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub.k27 %>%
         filter(batch == "New"), 
       aes(x = umap1, y = umap2, color = is.bad)) + 
  geom_point() + 
  facet_wrap(~ctype) + 
  theme_bw() + 
  geom_vline(xintercept = 15) +
  ggtitle("k27me3 new only") + 
  scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub.k27 %>%
         filter(batch == "Old"), 
       aes(x = umap1, y = umap2, color = is.bad)) + 
  geom_point() + 
  facet_wrap(~ctype) + 
  theme_bw() + 
  geom_vline(xintercept = 15) +
  ggtitle("k27me3 old") + 
  scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub.k27, aes(x = cell.var.within.sum.norm, fill = is.bad)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  ggtitle("k27me3 bad cells low var check") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# clean up NK cells 

ctype2col.k27 <- hash::hash(dat.meta.lst$k27me3$ctype.from.LL, dat.meta.lst$k27me3$colcode)

dat.meta.k27.clean <- jsub.k27 %>%
  filter(!is.bad)

ggplot(dat.meta.k27.clean, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k27me3 cleaned") + 
  scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                        guide = "legend") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.meta.k27.clean, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k27me3 cleaned") + 
  facet_wrap(~ctype.from.LL) + 
  scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                        guide = "legend") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

table(subset(dat.meta.k27.clean, batch == "New")$ctype)
table(subset(jsub.k27, batch == "New")$ctype)
table(subset(jsub.k27, batch == "New" & is.bad)$ctype)


# Clean up k9me3 ----------------------------------------------------------

ggplot(dat.meta.lst$k9me3, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k9me3") + 
  scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                        guide = "legend") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.lst$k9me3, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k9me3") + 
  scale_color_viridis_c(direction = -1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.lst$k9me3, aes(x = umap1, y = umap2, color = as.character(louvain))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k9me3") + 
  scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


jsub.k9me3 <- dat.meta.lst$k9me3 %>%
  mutate(is.bad = (louvain == "2" & umap1 > 12 & batch == "New"))

ggplot(jsub.k9me3,
       aes(x = umap1, y = umap2, color = is.bad)) + 
  geom_point() + 
  facet_wrap(~batch) + 
  theme_bw() + 
  geom_vline(xintercept = 12) +
  ggtitle("k9me3") + 
  scale_color_manual(values = cbPalette)  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub.k9me3, aes(x = cell.var.within.sum.norm, fill = is.bad)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  ggtitle("k9me3 bad cells low var check") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ctype2col.k9 <- hash::hash(dat.meta.lst$k9me3$ctype.from.LL, dat.meta.lst$k9me3$colcode)

dat.meta.k9.clean <- jsub.k9me3 %>%
  filter(!is.bad)

ggplot(dat.meta.k9.clean, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k9me3 cleaned") + 
  facet_wrap(~batch) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.k9.clean, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("k9me3 cleaned") + 
  facet_wrap(~ctype.from.LL) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

dat.meta.cleaned <- list(k27me3 = dat.meta.k27.clean,
                         k9me3 = dat.meta.k9.clean)


# Write count mats cleaned up  --------------------------------------------


for (jmark in jmarks){
  print(jmark)
  
  
  inf.dynbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed_dynamic_bins/BM_", jmark, "/count_mat_merged_with_old.", jmark, ".2022-01-26.rds")
  inf.allbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed/BM_", jmark, "/count_mat_merged_with_old.", jmark, ".2022-01-26.rds")
  assertthat::assert_that(file.exists(inf.dynbins))
  assertthat::assert_that(file.exists(inf.allbins))
  
  count.dynbins <- readRDS(inf.dynbins)
  count.allbins <- readRDS(inf.allbins)
  
  # filter bins
  jchromos.check <- sort(unique(sapply(rownames(count.allbins), function(x) strsplit(x, ":")[[1]][[1]])))
  print(jchromos.check)
  rnames.orig <- rownames(count.allbins)
  jchromos.auto <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  jchromos.grep <- paste(paste0("^", jchromos.auto), collapse = "|")
  rnames.i <- grepl(jchromos.grep, rnames.orig)
  rnames.filt <- rnames.orig[rnames.i]
  print(dim(count.allbins))
  count.allbins <- count.allbins[rnames.filt, ]
  print(dim(count.allbins))
  
  # load TSS
  if (jmark == "k27me3"){
    inf.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_rep2rep3only/count_mat_TSS_rep2rep3only.", jmark, ".2022-02-08.rds")
  } else {
    inf.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS/count_mat_TSS_combined.", jmark, ".2022-02-07.rds")
  }
  count.tss <- readRDS(inf.tss)
  
  
  metaout <- dat.meta.cleaned[[jmark]]
  
  cells.keep <- metaout$cell
  cells.keep.newonly <- subset(metaout, batch == "New")$cell
  
  count.dynbins <- readRDS(inf.dynbins)
  count.allbins <- readRDS(inf.allbins)
  count.tss <- readRDS(inf.tss)
  
  count.dynbins.filt <- count.dynbins[, cells.keep]
  count.allbins.filt <- count.allbins[, cells.keep]
  count.tss.filt <- count.tss[, cells.keep]
  
  count.dynbins.filt.newonly <- count.dynbins[, cells.keep.newonly]
  count.allbins.filt.newonly <- count.allbins[, cells.keep.newonly]
  count.tss.filt.newonly <- count.tss[, cells.keep.newonly]
  
  print(dim(count.dynbins))
  print(dim(count.dynbins.filt))
  
  print(dim(count.allbins))
  print(dim(count.allbins.filt))
  
  print(dim(count.tss))
  print(dim(count.tss.filt))
  
  outf.dynbins <- file.path(outdir, paste0("count_mat_cleaned_dynbins.", jmark, ".", Sys.Date(), ".rds"))
  outf.allbins <- file.path(outdir, paste0("count_mat_cleaned_allbins.", jmark, ".", Sys.Date(), ".rds"))
  outf.tss <- file.path(outdir, paste0("count_mat_cleaned_tss.", jmark, ".", Sys.Date(), ".rds"))
  
  outf.dynbins.newonly <- file.path(outdir, paste0("count_mat_cleaned_dynbins.newonly.", jmark, ".", Sys.Date(), ".rds"))
  outf.allbins.newonly <- file.path(outdir, paste0("count_mat_cleaned_allbins.newonly.", jmark, ".", Sys.Date(), ".rds"))
  outf.tss.newonly <- file.path(outdir, paste0("count_mat_cleaned_tss.newonly.", jmark, ".", Sys.Date(), ".rds"))
  
  outf.meta <- file.path(outdir, paste0("metadata_cleaned.", jmark, ".", Sys.Date(), ".txt"))
  
  # saveRDS(count.dynbins.filt, file = outf.dynbins)
  # saveRDS(count.allbins.filt, file = outf.allbins)
  # saveRDS(count.tss.filt, file = outf.tss)
  
  saveRDS(count.dynbins.filt.newonly, file = outf.dynbins.newonly)
  saveRDS(count.allbins.filt.newonly, file = outf.allbins.newonly)
  saveRDS(count.tss.filt.newonly, file = outf.tss.newonly)
  
  # fwrite(metaout, file = outf.meta)
  
}


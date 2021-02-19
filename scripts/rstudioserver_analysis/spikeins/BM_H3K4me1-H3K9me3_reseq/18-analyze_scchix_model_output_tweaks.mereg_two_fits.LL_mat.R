# Jake Yeung
# Date of Creation: 2021-02-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/18-analyze_scchix_model_output_tweaks.mereg_two_fits.LL_mat.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(ggforce)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/double_stain_outputs"
outpdf <- file.path(outdir, paste0("scchix_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".LL_mats.pdf"))
outrds <- file.path(outdir, paste0("scchix_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".act_repress_coord_merged_lst.rds"))

ctypes.order.k4 <- c("Eryths", "Bcells", "NKs", "DCs", "Granulocytes", "Basophils", "pDCs", "HSPCs")
names(ctypes.order.k4) <- ctypes.order.k4

ctypes.order.k9 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")
names(ctypes.order.k9) <- ctypes.order.k9

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


set.seed(0)
make.plots <- TRUE
write.tables <- FALSE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

# Load metas  -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})

# Load output -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

wmin <- 0.1
wmax <- 0.9
pcount <- 0

wmax.vec <- c(0.5, 0.9); names(wmax.vec) <- wmax.vec

# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.4.wmax_0.6.RData")
# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.wmax_0.7.RData")
# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.wmax_0.7.pseudocount_0.RData")

infsrdata <- lapply(wmax.vec, function(wmax){
  print(wmax)
  infrdata <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_", wmin, ".wmax_", wmax, ".pseudocount_", pcount, ".RData"))
  assertthat::assert_that(file.exists(infrdata))
  return(infrdata)
})

fits.out.lst <- lapply(infsrdata, function(infrdata){
  load(infrdata, v=T)
  fits.out <- act.repress.coord.lst
  w.lst <- sapply(fits.out, function(x) x$w)
  
  # if louvains are now from clusters need eto rethink jcoord
  cell.vec <- names(fits.out)
  names(cell.vec) <- cell.vec
  coords.dbl <- lapply(cell.vec, function(jcell){
    jfit <- fits.out[[jcell]]
    jweight <- fits.out[[jcell]]$w
    p.mat <- SoftMax(jfit$ll.mat)
    jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
    jmax <- max(p.mat)
    
    # rows are active, columns are repress I THINK?
    # TODO: assumes underscores be careful!
    jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
    jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
    
    if (grepl("_", jlouv.act)){
      jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
    }
    if (grepl("_", jlouv.repress)){
      jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
    }
    out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
    return(out.dat)
  }) %>%
    bind_rows()
  return(coords.dbl)
})

act.repress.coord.lst.lst <- lapply(infsrdata, function(infrdata){
  load(infrdata, v=T)
  return(act.repress.coord.lst)
})



# Check erythroblast estimates for both models  ---------------------------

jctype <- "Basophils"
jctype <- "pDCs"
jctype <- "HSPCs"
jctype <- "Eryths"

jctypes <- unique(fits.out.lst[[1]]$louv.act)
for (jctype in jctypes){
  dat.ctype1 <- subset(fits.out.lst[[1]], louv.act == jctype)
  dat.ctype2 <- subset(fits.out.lst[[2]], louv.act == jctype)
  print(jctype)
  print(table(dat.ctype1$louv.repress))
  print(table(dat.ctype2$louv.repress))
}


# Mix in fits  ------------------------------------------------------------

fits.eryth <- subset(fits.out.lst[["0.5"]], louv.act == "Eryths") %>%
  mutate(wmax = 0.5)
fits.noneryth <- subset(fits.out.lst[["0.9"]], louv.act != "Eryths") %>%
  mutate(wmax = 0.9)

fits.merge <- rbind(fits.noneryth, fits.eryth)

cells.eryths <- subset(fits.merge, louv.act == "Eryths")$cell
cells.noeryths <- subset(fits.merge, louv.act != "Eryths")$cell

act.repress.coord.eryths.lst <- act.repress.coord.lst.lst$`0.5`[cells.eryths]
act.repress.coord.noeryths.lst <- act.repress.coord.lst.lst$`0.9`[cells.noeryths]

act.repress.coord.merge.lst <- c(act.repress.coord.eryths.lst, act.repress.coord.noeryths.lst)

# Plot some LLs  ----------------------------------------------------------

jclsts <- unique(fits.merge$louv.act); names(jclsts) <- jclsts


# pdf(outpdf, useDingbats = FALSE)
for (jclst in jclsts){
  print(jclst)
  if (jclst == "Basophils"){
    jsub.tmp <- subset(fits.merge, louv.act == jclst & louv.repress == "Granulocytes" & w < (wmax - 0.01)) %>%
      arrange(lnprob)
  } else {
    jsub.tmp <- subset(fits.merge, louv.act == jclst & lnprob == 0 & w < (wmax - 0.01))
  }
  jcell.tmp <- jsub.tmp$cell[[1]]
  
  LL.mat.tmp <- act.repress.coord.merge.lst[[jcell.tmp]]$ll.mat
  LL.mat.tmp <- data.frame(K4me1 = rownames(LL.mat.tmp), LL.mat.tmp, stringsAsFactors = FALSE)
  
  LL.long.tmp <- data.table::melt(LL.mat.tmp, id.vars = "K4me1", value.name = "LL", variable.name = "K9me3")
  LL.long.tmp$K4me1 <- factor(LL.long.tmp$K4me1, levels = ctypes.order.k4)
  LL.long.tmp$K9me3 <- factor(LL.long.tmp$K9me3, levels = ctypes.order.k9)
  
  m <- ggplot(LL.long.tmp, aes(x = K4me1, y = K9me3, fill = LL)) + 
    geom_tile() + 
    theme_bw() + 
    scale_fill_viridis_c() + 
    ggtitle(paste(jclst, "cell", jcell.tmp)) + 
    theme(aspect.ratio= 4/8, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          legend.position = "bottom")
  print(m)
}
dev.off()

saveRDS(act.repress.coord.merge.lst, file = outrds)


# ggplot(LL.long.tmp, aes(x = K4me1, y = K9me3, fill = -LL)) + 
#   geom_tile() + 
#   theme_bw() + 
#   scale_fill_viridis_c() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


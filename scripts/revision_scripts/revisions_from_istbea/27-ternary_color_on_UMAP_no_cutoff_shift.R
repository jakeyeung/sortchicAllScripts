# Jake Yeung
# Date of Creation: 2022-04-26
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/27-ternary_color_on_UMAP.R
# Add three colors to UMAP 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(tricolore)
library(parallel)

ncores <- 4

EuclideanDistance <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs"
outpdf <- file.path(outdir, paste0("imputation_sca1_kit_lin_no_cutoff_shift_sl4_only.", Sys.Date(), ".pdf"))
outrdata <- file.path(outdir, paste0("imputation_sca1_kit_lin_no_cutoff_shift_sl4_only.", Sys.Date(), ".RData"))
pdf(outpdf, useDingbats = FALSE)

# Load PCA for euclidean distance  --------------------------------------------------------------


indir.impute <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects"

dat.impute.lst <- lapply(jmarks, function(jmark){
  inf.impute <- file.path(indir.impute, paste0("dat_impute_bins_", jmark, ".2022-04-24.rds"))
  dat.impute <- readRDS(inf.impute)
})

# # indir.pca <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"
# dat.pca.lst <- lapply(jmarks, function(jmark){
#   inf.pca <- file.path(indir.pca, paste0("pca_output.", jmark, ".rds"))
#   dat.pca <- readRDS(inf.pca)
# })



# Load meta ----------------------------------------------------------------


experis.keep <- paste(c("PZ-sortChIC-BM-SL3", "PZ-sortChIC-BM-SL4", "PZ-sortChIC-BM-SL5"), collapse = "|")

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    rowwise() %>%
    mutate(experi = ClipLast(ClipLast(cell, jsep = "_"), jsep = "-"),
           experi.keep = grepl(experis.keep, experi))
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
dat.meta.colors.add <- data.frame(ctype.from.LL = c("GMP", "Tcells"), colcode = c("#d88543", "#04A804"))
dat.meta.colors <- rbind(dat.meta.colors, dat.meta.colors.add)

ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)



# Load facs info  ---------------------------------------------------------

dat.facs.lst <- lapply(jmarks, function(jmark){
  inf.facs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/from_peter/metadata_inclFACS.", jmark, ".2022-03-01.txt")
  jdat <- as.data.frame(fread(inf.facs))
  colnames(jdat) <- gsub("-", "_", colnames(jdat))
  jgrep <- paste(c("kit", "^Sca", "^lin"), collapse = "|")
  cnames.keep <- c("cell", grep(jgrep, colnames(jdat), value = TRUE))
  jdat <- jdat[, cnames.keep]
  jdat <- jdat[, sort(colnames(jdat))]
  return(jdat)
})

print(head(dat.facs.lst[[1]]))

# Add to UMAP  ------------------------------------------------------------

dat.annot.lst <- lapply(jmarks, function(jmark){
  left_join(dat.meta.lst[[jmark]], dat.facs.lst[[jmark]])
})              

jcheck.dup.lst <- lapply(jmarks, function(jmark){
  jcheck <- dat.annot.lst[[jmark]] %>% filter(!is.na(Sca_1_PeCy7_log2)) 
  jcheck <- dat.annot.lst[[jmark]] %>% 
    rowwise() %>%
    mutate(ckit.na = is.na(C_kit_BB700), 
           lin.na = is.na(lin_PE), 
           sca1.na = is.na(Sca_1_PeCy7)) %>%
    dplyr::select(experi, batch, ckit.na, lin.na, sca1.na) %>%
    filter(batch == "New") %>%
    filter(!(ckit.na & lin.na & sca1.na)) 
  
  jcheck.dup <- jcheck[!duplicated(jcheck), ]
})

jcheck.dup.lst



# Average out Sca1 and Ckit -----------------------------------------------


dat.sub.lst <- lapply(dat.annot.lst, function(jdat){
  subset(jdat, grepl(experis.keep, cell))
})


# Impute missing values ---------------------------------------------------

dat.sub.impute.lst <- lapply(dat.sub.lst, function(jdat){
  jdat <- jdat %>%
    ungroup() %>%
    mutate(ckit_zscore = scale(C_kit_BB700_log2), 
           sca1_zscore = scale(Sca_1_PeCy7_log2),
           lin_zscore = scale(lin_PE_log2),
           ckit_impute = ifelse(is.na(ckit_zscore), min(ckit_zscore, na.rm = TRUE), ckit_zscore), 
           sca1_impute = ifelse(is.na(sca1_zscore), min(sca1_zscore, na.rm = TRUE), sca1_zscore),
           lin_impute = ifelse(is.na(lin_zscore), min(lin_zscore, na.rm = TRUE), lin_zscore)) %>%
    rowwise() %>%
    mutate(ckit_f = SoftMax(x = c(ckit_impute, sca1_impute, lin_impute), return.log = FALSE, logfn = log)[[1]],
           sca1_f = SoftMax(x = c(ckit_impute, sca1_impute, lin_impute), return.log = FALSE, logfn = log)[[2]],
           lin_f = SoftMax(x = c(ckit_impute, sca1_impute, lin_impute), return.log = FALSE, logfn = log)[[3]]) %>%
    filter( !(is.na(ckit_zscore) & is.na(sca1_zscore) & is.na(lin_zscore)) )
})


data.frame(subset(dat.annot.lst$k4me1, cell == "PZ-sortChIC-BM-SL4-k4me1-2_90"))

jmark <- "k4me3"
jmark <- "k4me1"

for (jmark in jmarks){
  
  m <- ggplot(dat.annot.lst[[jmark]] %>% arrange(desc(is.na(Sca_1_PeCy7_log2))) %>%
           filter(batch == "New"), 
         aes(x = umap1, y = umap2, color = Sca_1_PeCy7_log2)) + 
    geom_point() + 
    scale_color_viridis_c(na.value = "grey85") + 
    ggtitle(jmark) + 
    theme_bw() + 
    facet_wrap(~experi) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.annot.lst[[jmark]] %>% arrange(desc(is.na(Sca_1_PeCy7_log2))),
         aes(x = umap1, y = umap2, color = Sca_1_PeCy7_log2)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_viridis_c(na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.annot.lst[[jmark]] %>% arrange(desc(is.na(C_kit_BB700_log2))) %>%
           filter(batch == "New"), 
         aes(x = umap1, y = umap2, color = C_kit_BB700_log2)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_viridis_c(na.value = "grey85") + 
    facet_wrap(~experi) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.annot.lst[[jmark]] %>% arrange(desc(is.na(C_kit_BB700_log2))),
         aes(x = umap1, y = umap2, color = C_kit_BB700_log2)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_viridis_c(na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.annot.lst[[jmark]] %>% arrange(desc(is.na(lin_PE_log2))) %>%
           filter(batch == "New"), 
         aes(x = umap1, y = umap2, color = lin_PE_log2)) + 
    geom_point() + 
    ggtitle(jmark) + 
    facet_wrap(~experi) + 
    scale_color_viridis_c(na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.annot.lst[[jmark]] %>% arrange(desc(is.na(lin_PE_log2))),
         aes(x = umap1, y = umap2, color = lin_PE_log2)) + 
    geom_point() + 
    ggtitle(jmark) + 
    facet_wrap(~experi) + 
    scale_color_viridis_c(na.value = "grey85") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
}



# Get distances between cells  --------------------------------------------

print(jcheck.dup.lst[[1]])

cell2label.lin.lst <- lapply(dat.sub.impute.lst, function(jdat){
  hash::hash(jdat$cell, jdat$lin_zscore)
})

cell2label.sca1.lst <- lapply(dat.sub.impute.lst, function(jdat){
  hash::hash(jdat$cell, jdat$sca1_zscore)
})





# Use sl4 only  -----------------------------------------------------------


mat.sl4.lst <- lapply(jmarks, function(jmark){
  cells.keep <- dat.sub.impute.lst[[jmark]]$cell
  mat <- dat.impute.lst[[jmark]]
  cnames.keep.i <- grepl("SL4", colnames(mat))
  cnames.keep <- colnames(mat)[cnames.keep.i]
  cells.keep2 <- cnames.keep[cnames.keep %in% cells.keep]
  return(mat[, cells.keep2])
})



# sl3 needs to impute from sl4/sl5 (sl4/sl5 has lin)


mat.sl4sl5.lst <- lapply(jmarks, function(jmark){
  cells.keep <- dat.sub.impute.lst[[jmark]]$cell
  mat <- dat.impute.lst[[jmark]]
  cnames.keep.i <- grepl("SL4|SL5", colnames(mat))
  cnames.keep <- colnames(mat)[cnames.keep.i]
  cells.keep2 <- cnames.keep[cnames.keep %in% cells.keep]
  return(mat[, cells.keep2])
})

# sl5 needs to impute from sl3/sl4 (sl3/sl4 has sca1)
mat.sl3sl4.lst <- lapply(jmarks, function(jmark){
  cells.keep <- dat.sub.impute.lst[[jmark]]$cell
  mat <- dat.impute.lst[[jmark]]
  cnames.keep.i <- grepl("SL3|SL4", colnames(mat))
  cnames.keep <- colnames(mat)[cnames.keep.i]
  cells.keep2 <- cnames.keep[cnames.keep %in% cells.keep]
  return(mat[, cells.keep2])
})

# check there's a label for each cell 
sl4sl5.check.lst <- lapply(jmarks, function(jmark){
  sapply(colnames(mat.sl4sl5.lst[[jmark]]), AssignHash, cell2label.lin.lst[[jmark]], null.fill = NA)
})
lapply(sl4sl5.check.lst, function(jlst) range(jlst))

sl3sl4.check.lst <- lapply(jmarks, function(jmark){
  sapply(colnames(mat.sl3sl4.lst[[jmark]]), AssignHash, cell2label.sca1.lst[[jmark]], null.fill = NA)
})
lapply(sl3sl4.check.lst, function(jlst) range(jlst))

# Get nearest for SL3 and SL5 ---------------------------------------------

cells.sl3.lst <- lapply(dat.impute.lst, function(mat){
  cnames.keep <- grepl("SL3", colnames(mat))
  cells <- colnames(mat)[cnames.keep]
  names(cells) <- cells
  return(cells)
})

cells.sl5.lst <- lapply(dat.impute.lst, function(mat){
  cnames.keep <- grepl("SL5", colnames(mat))
  cells <- colnames(mat)[cnames.keep]
  names(cells) <- cells
  return(cells)
})

# SL3 need to impute lineage 
keepn <- 10

library(parallel)

# bad cell PZ-sortChIC-BM-SL3-k4me1-2_140

cells.sl3.nearest.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jcells <- cells.sl3.lst[[jmark]]
  closest.cells.labels.lst <- parallel::mclapply(jcells, function(jcell){
    
    ref.vec <- dat.impute.lst[[jmark]][, jcell]
    # euc.dists <- apply(mat.sl4sl5.lst[[jmark]], 2, function(jcol){
    euc.dists <- apply(mat.sl4.lst[[jmark]], 2, function(jcol){
      EuclideanDistance(jcol, ref.vec)
    }) %>%
      sort(decreasing = FALSE)
    # get labels for top N
    closest.cells <- names(euc.dists[1:keepn])
    closest.cells.labels <- sapply(closest.cells, function(x){
      AssignHash(x = x, jhash = cell2label.lin.lst[[jmark]], null.fill = NA)
    })
    assertthat::assert_that(all(!is.na(closest.cells.labels)))
    return(closest.cells.labels)
  }, mc.cores = ncores)
  return(closest.cells.labels.lst)
})

# sl5 cells need sca1

# bad cell  
# jcell.check <- "PZ-sortChIC-BM-SL5-k27me3-3_15"
cells.sl5.nearest.lst.out <- lapply(jmarks, function(jmark){
  print(jmark)
  jcells <- cells.sl5.lst[[jmark]]
  closest.cells.labels.lst <- parallel::mclapply(jcells, function(jcell){
    ref.vec <- dat.impute.lst[[jmark]][, jcell]
    # euc.dists <- apply(mat.sl3sl4.lst[[jmark]], 2, function(jcol){
    euc.dists <- apply(mat.sl4.lst[[jmark]], 2, function(jcol){
      EuclideanDistance(jcol, ref.vec)
    }) %>%
      sort(decreasing = FALSE)
    # get labels for top N
    closest.cells <- names(euc.dists[1:keepn])
    closest.cells.labels <- sapply(closest.cells, function(x){
      AssignHash(x = x, jhash = cell2label.sca1.lst[[jmark]], null.fill = NA)
    })
    assertthat::assert_that(all(!is.na(closest.cells.labels)))
    return(list(closest.cells.labels = closest.cells.labels, euc.dists = euc.dists[1:keepn]))
  }, mc.cores = ncores)
  return(closest.cells.labels.lst)
})

cells.sl5.nearest.dists.lst <- lapply(cells.sl5.nearest.lst.out, function(jcells){
  lapply(jcells, function(jcell){
    jcell$euc.dists
  })
})

cells.sl5.nearest.lst <- lapply(cells.sl5.nearest.lst.out, function(jcells){
  lapply(jcells, function(jcell){
    jcell$closest.cells.labels
  })
})

# plot closest
# cells.sl5.nearest.dists1.lst <- cells.sl5.nearest.dists.lst
plot(density(sapply(cells.sl5.nearest.dists.lst$k4me1, function(x) x[[1]])))
plot(density(sapply(cells.sl5.nearest.dists.lst$k4me3, function(x) x[[1]])))
plot(density(sapply(cells.sl5.nearest.dists.lst$k27me3, function(x) x[[1]])))
plot(density(sapply(cells.sl5.nearest.dists.lst$k9me3, function(x) x[[1]])))

dist.means.lst <- lapply(cells.sl5.nearest.dists.lst, function(eucdists){
  jmean <- mean(sapply(eucdists, function(x) x[[1]]))
  jsd <- sd(sapply(eucdists, function(x) x[[1]]))
  jmean + 0 * jsd
})

plot(dat.sub.impute.lst$k4me1$lin_zscore, dat.sub.impute.lst$k4me1$sca1_zscore)
plot(dat.sub.impute.lst$k4me1$lin_zscore, dat.sub.impute.lst$k4me1$ckit_zscore)
plot(dat.sub.impute.lst$k4me1$sca1_zscore, dat.sub.impute.lst$k4me1$ckit_zscore)


dat.dist.lst <- lapply(jmarks, function(jmark){
  data.frame(cell = names(cells.sl5.nearest.dists.lst[[jmark]]), 
             eucdist = sapply(cells.sl5.nearest.dists.lst[[jmark]], function(x) x[[1]]), stringsAsFactors = FALSE)
})

dat.dist.annot.lst <- lapply(jmarks, function(jmark){
  left_join(dat.dist.lst[[jmark]], dat.sub.impute.lst[[jmark]]) %>%
    rowwise() %>%
    mutate(beyond.cutoff = eucdist > dist.means.lst[[jmark]])
})

dat.beyond.cutoff.lst <- lapply(dat.dist.annot.lst, function(jdat){
  data.frame(cell = jdat$cell, beyond.cutoff = jdat$beyond.cutoff, stringsAsFactors = FALSE)
})

# plot on UMAP 
m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.dist.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = beyond.cutoff)) + 
    geom_point() + 
    # scale_color_viridis_c() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})
print(m.lst)

# Mean and variance of closest cells  -------------------------------------

cells.sl3.nearest.mean.lst <- lapply(cells.sl3.nearest.lst, function(jcells.lst){
  lapply(jcells.lst, function(jcells){
    mean(jcells)
  })
})

cells.sl3.nearest.var.lst <- lapply(cells.sl3.nearest.lst, function(jcells.lst){
  lapply(jcells.lst, function(jcells){
    var(jcells)
  })
})


cells.sl5.nearest.mean.lst <- lapply(cells.sl5.nearest.lst, function(jcells.lst){
  lapply(jcells.lst, function(jcells){
    mean(jcells)
  })
})

cells.sl5.nearest.var.lst <- lapply(cells.sl5.nearest.lst, function(jcells.lst){
  lapply(jcells.lst, function(jcells){
    var(jcells)
  })
})

# Impute KNN  -------------------------------------------------------------

save.image(file = outrdata)

dat.sub.impute.knn.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.sub.impute.lst[[jmark]] %>%
    left_join(., dat.beyond.cutoff.lst[[jmark]]) %>%
    mutate(beyond.cutoff = ifelse(is.na(beyond.cutoff), FALSE, beyond.cutoff))
  # ckit does not need imputing just assign zscore
  sca1_zscore_min <- min(jdat$sca1_zscore, na.rm = TRUE)
  jdat <- jdat %>%
    rowwise() %>%
    mutate(lin_knn_impute = ifelse(is.na(lin_zscore), cells.sl3.nearest.mean.lst[[jmark]][[cell]], lin_zscore), 
           sca1_knn_impute = ifelse(is.na(sca1_zscore), cells.sl5.nearest.mean.lst[[jmark]][[cell]], sca1_zscore), 
           ckit_knn_impute = ckit_zscore) %>%
    rowwise() %>%
    # mutate(sca1_knn_impute = ifelse(beyond.cutoff, sca1_zscore_min, sca1_knn_impute)) %>%
    mutate(sca1_knn_impute = sca1_knn_impute - 0.5) %>%
    mutate(ckit_f_impute = SoftMax(x = c(ckit_knn_impute, sca1_knn_impute, lin_knn_impute), return.log = FALSE, logfn = log)[[1]],
           sca1_f_impute = SoftMax(x = c(ckit_knn_impute, sca1_knn_impute, lin_knn_impute), return.log = FALSE, logfn = log)[[2]],
           lin_f_impute = SoftMax(x = c(ckit_knn_impute, sca1_knn_impute, lin_knn_impute), return.log = FALSE, logfn = log)[[3]]) 
  # range(jdat$lin_knn_impute)
  return(jdat)
})


# Downstream impute -------------------------------------------------------


jmark <- "k4me1"
ggplot(dat.sub.impute.knn.lst[[jmark]] %>% arrange(desc(is.na(sca1_knn_impute))) %>%
         filter(batch == "New"), 
       aes(x = umap1, y = umap2, color = sca1_f_impute)) + 
  geom_point() + 
  scale_color_viridis_c(na.value = "grey85") + 
  theme_bw() + 
  facet_wrap(~experi) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.sub.impute.knn.lst[[jmark]] %>% arrange(desc(is.na(ckit_impute))) %>%
         filter(batch == "New"), 
       aes(x = umap1, y = umap2, color = ckit_f_impute)) + 
  geom_point() + 
  scale_color_viridis_c(na.value = "grey85") + 
  facet_wrap(~experi) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.sub.impute.knn.lst[[jmark]] %>% arrange(desc(is.na(lin_impute))) %>%
         filter(batch == "New"), 
       aes(x = umap1, y = umap2, color = lin_f_impute)) + 
  geom_point() + 
  facet_wrap(~experi) + 
  scale_color_viridis_c(na.value = "grey85") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Ternary  ----------------------------------------------------------------

# jmark <- "k9me3"
jmark <- "k4me3"

jmark <- "k9me3"
jmark <- "k27me3"
jmark <- "k4me3"

jmark <- "k4me1"
for (jmark in jmarks){
  
  jsub <- dat.sub.impute.knn.lst[[jmark]]
  colors_and_legend <- tricolore::Tricolore(jsub, 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
  jsub$rgb <- colors_and_legend$rgb
  
  print(colors_and_legend$key)
  # m <- ggplot(jsub, aes(x = umap1, y = umap2, color = rgb)) + 
  #   geom_point() + 
  #   scale_color_identity() + 
  #   ggtitle(jmark) + 
  #   theme_bw() + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # add other cells as grey underneat
  m <- ggplot() + 
    geom_point(mapping = aes(x = umap1, umap2), data = dat.meta.lst[[jmark]], color = "grey85") + 
    geom_point(mapping = aes(x = umap1, umap2, color = rgb), data = jsub) + 
    scale_color_identity() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  print(m)
}


# Do all  -----------------------------------------------------------------


# pdf("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/imputation_outputs/all_triangle.pdf", useDingbats = FALSE)

jsub.all <- lapply(jmarks, function(jmark){ 
  dat.sub.impute.knn.lst[[jmark]] %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows() %>%
  dplyr::rename(Sca1 = "sca1_f_impute",
                cKit = "ckit_f_impute", 
                Lin = "lin_f_impute")

colors_and_legend <- tricolore::Tricolore(jsub.all, 'Sca1', 'cKit', 'Lin')
jsub.all$rgb <- colors_and_legend$rgb
print(colors_and_legend$key)

# dev.off()



dev.off()

save.image(file = outrdata)


# 
# # Impute using norm -------------------------------------------------------
# 
# library(norm2)
# 
# jmark <- "k4me1"
# mat.input <- dat.sub.impute.lst[[jmark]] %>%
#   dplyr::select(cell, c(ckit_zscore, sca1_zscore, lin_zscore)) %>%
#   as.data.frame()
# 
# rownames(mat.input) <- mat.input$cell
# mat.input$cell <- NULL
# 
# em.out <- norm2::emNorm(mat.input, prior="ridge", prior.df=1/20)
# summary(em.out)
# 
# impute.out <- impNorm(em.out)
# 
# colnames(impute.out) <- paste(colnames(impute.out), "impute", sep = "_")
# 
# # colnames
# 
# # show on UMAP 
# dat.impute <- data.frame(cell = rownames(impute.out), impute.out, stringsAsFactors = FALSE) %>%
#   rowwise() %>%
#     mutate(ckit_f_impute = SoftMax(x = c(ckit_zscore_impute, sca1_zscore_impute, lin_zscore_impute), return.log = FALSE, logfn = log)[[1]],
#            sca1_f_impute = SoftMax(x = c(ckit_zscore_impute, sca1_zscore_impute, lin_zscore_impute), return.log = FALSE, logfn = log)[[2]],
#            lin_f_impute = SoftMax(x = c(ckit_zscore_impute, sca1_zscore_impute, lin_zscore_impute), return.log = FALSE, logfn = log)[[3]])  %>%
#   left_join(., dat.sub.impute.lst[[jmark]][, c("cell", "umap1", "umap2")])
# 
# ggplot(dat.impute, aes(x = umap1, y = umap2, color = sca1_zscore_impute)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_viridis_c() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(dat.impute, aes(x = umap1, y = umap2, color = sca1_f_impute)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_viridis_c() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# 
# ggplot(dat.sub.impute.knn.lst[[jmark]] %>% arrange(desc(is.na(sca1_knn_impute))) %>%
#          filter(batch == "New"), 
#        aes(x = umap1, y = umap2, color = sca1_f_impute)) + 
#   geom_point() + 
#   scale_color_viridis_c(na.value = "grey85") + 
#   theme_bw() + 
#   facet_wrap(~experi) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# ggplot(dat.sub.impute.knn.lst[[jmark]] %>% arrange(desc(is.na(sca1_knn_impute))) %>%
#          filter(batch == "New"), 
#        aes(x = umap1, y = umap2, color = sca1_knn_impute)) + 
#   geom_point() + 
#   scale_color_viridis_c(na.value = "grey85") + 
#   theme_bw() + 
#   facet_wrap(~experi) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 

# 
# # Test on Sl4 -------------------------------------------------------------
# 
# 
# jdat.sub.impute.filt.lst <- lapply(dat.sub.impute.lst, function(jdat){
#   as.data.frame(subset(jdat, grepl("SL4", cell), select = c(cell, umap1, umap2, sca1_impute, lin_impute, ckit_impute)))
# })
# 
# ggplot(jdat.sub.impute.filt.lst[[jmark]] %>% arrange(desc(is.na(lin_impute))),
#        aes(x = umap1, y = umap2, color = lin_impute)) + 
#   geom_point() + 
#   scale_color_viridis_c(na.value = "grey85") + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# library(tricolore)
# 
# jsub <- jdat.sub.impute.filt.lst[[1]]
# P <- as.data.frame(prop.table(matrix(runif(3^6), ncol = 3), 1))
# colors_and_legend <- tricolore::Tricolore(jsub, 'sca1_impute', 'ckit_impute', 'lin_impute', p1 = 'ed_0to2', p2 = 'ed_3to4', p3 = 'ed_5to8')
# head(colors_and_legend$rgb)
# 
# jsub$rgb <- colors_and_legend$rgb
# 
# colors_and_legend$key
# 
# ggplot(jsub, aes(x = umap1, y = umap2, color = rgb)) + 
#   geom_point() + 
#   scale_color_identity() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # DemoTricolore()

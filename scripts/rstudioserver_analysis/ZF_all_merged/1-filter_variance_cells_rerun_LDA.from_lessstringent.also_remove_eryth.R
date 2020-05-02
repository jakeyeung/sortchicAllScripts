# Jake Yeung
# Date of Creation: 2020-04-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/1-filter_variance_cells_rerun_LDA.R
# Filte variance rerun LDA 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(gtools)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmark <- "H3K4me3"
jmark <- "H3K4me1"
jmark <- "H3K9me3"
jmark <- "H3K27me3"

jwin <- 50000L
jprefix <- "/home/jyeung/hub_oudenaarden/jyeung"
jsuffix <- ".lessstringent"
jsuffix2 <- "remove_eryth"
indir.lda <- file.path(jprefix, paste0("data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_", jwin, jsuffix))


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jcutoffs <- c(0.75, 2, 1, 0.5); names(jcutoffs) <- jmarks

bad.clsts <- list(H3K4me1 = c("10", "6"), 
                  H3K4me3 = c("2"),
                  H3K27me3 = c("3", "5", "9"),
                  H3K9me3 = c("1", "4"))

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/check_var_filt_cells"
outdir <- file.path(jprefix, paste0("data/zebrafish_scchic/from_rstudio/impute_var_filt_from", jsuffix, ".", jsuffix2))
dir.create(outdir)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


for (jmark in jmarks){
  outpdf <- file.path(outdir, paste0("plots.", jmark, ".umap_var_checks.pdf"))
  outrds <- file.path(outdir, paste0("counts_table_var_filt.", jmark, ".imputevar_", jcutoffs[[jmark]], ".rds"))
  
  # Load TSS signal ---------------------------------------------------------
  
  fname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
  inf.lda <- file.path(indir.lda, fname)
  print(inf.lda)
  
  # Load LDA ----------------------------------------------------------------
  
  load(inf.lda, v=T)
  
  rows.reorder <- mixedorder(rownames(count.mat))
  
  count.mat.sorted <- count.mat[rows.reorder, ]
  
  # Do UMAP  ----------------------------------------------------------------
  
  jchromos <- paste("chr", seq(25), sep = "")
  
  tm.result <- posterior(out.lda)
  topics.mat <- tm.result$topics
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.umap <- DoUmapAndLouvain(topics.mat, jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"))
  
  m.before <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point() +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark) + scale_color_manual(values = cbPalette)
  
  # filter bad clusters
  bad.clst <- bad.clsts[[jmark]]
  print("Clusters before:")
  print(unique(dat.umap$louvain))
  dat.umap <- subset(dat.umap, !louvain %in% bad.clst)
  print("Clusters after:")
  print(unique(dat.umap$louvain))
  
  dat.var <- CalculateVarAll(dat.impute.log[rows.reorder, ], jchromos)
  
  dat.var.raw <- CalculateVarRaw(count.mat.sorted, merge.size = 100, jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)
  
  dat.merge <- left_join(dat.umap, dat.var) %>%
    left_join(., dat.var.raw)
  
  m1 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() +  
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
 
  m2 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = log10(ncuts.var))) + geom_point() +  
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  
  m3 <- ggplot(dat.merge, aes(x = ncuts.var, y = cell.var.within.sum.norm, color = log10(ncuts))) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + scale_color_viridis_c() + scale_x_log10() + scale_y_log10() 
  
  # Filter by imputed variance ----------------------------------------------
  
  dat.merge <- dat.merge %>%
    rowwise() %>%
    mutate(low.var = cell.var.within.sum.norm < jcutoffs[[jmark]])
  
  m4 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = low.var)) + geom_point() +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark)
  
  m5 <- ggplot(dat.merge, aes(x = cell.var.within.sum.norm, fill = plate)) + geom_density(alpha = 0.25)  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  m6 <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + geom_point() +  
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~plate) + ggtitle(jmark) + scale_color_manual(values = cbPalette)
  
  pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
    print(m1)
    print(m2)
    print(m3)
    print(m4)
    print(m5)
    print(m6)
    print(m.before)
  dev.off()
  
  # write filtered mat to output
  cells.keep <- subset(dat.merge, !low.var)$cell
  count.mat.filt <- count.mat.sorted[, cells.keep]
  saveRDS(count.mat.filt, file = outrds)
}



# Jake Yeung
# Date of Creation: 2021-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/FACS/explore_FACS.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


WellToCellNumber <- function(w, totalcols = 24, shift = -1){
  # Well name (e.g., P20 to Plate coordinate (e.g. 380)
  jlet <- toupper(letters[1:26])
  jrow <- as.numeric(match(gsub("[^a-zA-Z]", "", w), jlet))
  jcol <- as.numeric(gsub("[a-zA-Z]", "", w))
  # shift -1 gives zerobased coordinates
  return( as.character((jrow - 1) * totalcols + jcol + shift) )
}



# Load metas  -------------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  fread(inf.meta)
})


# Load FACS  --------------------------------------------------------------


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
jplate <- "PZ-BM-rep3-H3K27me3-3"
jplates <- unique(subset(dat.metas$H3K27me3, jrep == "rep3")$experi)
names(jplates) <- jplates

cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME", "Tray X", "Tray Y", "Trigger Pulse Width",
                   "Time 1", "Time 2", "Time 3", "Drop Phase", "ROI Bits 1-16", "ROI Bits 17-32", "Classifier Bits", "Sort Result Bits",
                   "Sort Enable Bits")

# "*[561] 585/29" "*[640] 670/30" "*[561] 750 LP"

cnames.switchlist <- list("X..561..585.29" = "Lineage",
                          "X..640..670.30" = "Ckit",
                          "X..561..750.LP" = "Sca1",
                          "*[561] 585/29" = "Lineage",
                          "*[561] 610/20" = "Lineage",
                          "*[488] 710/50" = "Lineage",
                          "*[640] 670/30" = "Ckit",
                          "*[561] 750 LP" = "Sca1")


# cname.lineage <- "*[488] 710/50"
cname.lineage <- "*[561] 585/29"

# cname.lineage <- "*[561] 610/20"
cname.ckit <- "*[640] 670/30"
cname.lsk <- "*[561] 750 LP"
cnames.sorting <- c(cname.lineage, cname.ckit, cname.lsk)

# names(cnames.switchlist) <- cnames.switchlist
cnames.rename <- hash::hash(cnames.switchlist)
  
jmerge.lst <- lapply(jplates, function(jplate){
  print(jplate)
  
  inf.test <- file.path(hubprefix, paste0("jyeung/data/scChiC/facs_index_data/BM_round_3/H3K4me1_H3K27me3_20200901/", jplate, "_index.csv"))
  assertthat::assert_that(file.exists(inf.test))
  
  dat.facs <- fread(inf.test)
  
  dat.facs <- dat.facs[!duplicated(dat.facs$Well), ]
  
  # remove weird colnames
  
  
  # cnames.keep <- !colnames(dat.facs) %in% cnames.remove
  cnames.keep <- cnames.sorting
  print(cnames.keep)
  dat.facs.mat <- as.data.frame(dat.facs[, ..cnames.keep])
  
  cellnumbr <- sapply(dat.facs$Well, function(x) WellToCellNumber(x, totalcols = 24, shift = -1))
  rownames(dat.facs.mat) <- cellnumbr
  feats.keep <- which(apply(dat.facs.mat, 2, var) > 0)
  # dat.facs.mat.filt <- scale(log2(dat.facs.mat[, feats.keep]), center = TRUE, scale = TRUE)
  dat.facs.mat.filt <- dat.facs.mat[, feats.keep]
  
  for (cname.tmp in colnames(dat.facs.mat.filt)){
    cname.i <- which(colnames(dat.facs.mat.filt) == cname.tmp)
    cname.new <- cnames.switchlist[[cname.tmp]]
    colnames(dat.facs.mat.filt)[cname.i] <- cname.new
  }
  cnames.new <- colnames(dat.facs.mat.filt)
  
  
  # dat.facs.umap <- DoUmapAndLouvain(topics.mat = dat.facs.mat.filt, jsettings = jsettings)
  dat.facs.umap <- as.data.frame(dat.facs.mat.filt)
  dat.facs.umap$cellplate <- paste(jplate, rownames(dat.facs.mat), sep = "_")
  dat.facs.umap$platename <- jplate
  
  # Merge with known cell types  --------------------------------------------
  
  jmerge <- left_join(dat.metas$H3K27me3, dat.facs.umap, by = c("cell" = "cellplate"))
  # jmerge2 <- left_join(jmerge, data.frame(cellplate = paste(jplate, rownames(dat.facs.mat.filt), sep = "_"), dat.facs.mat.filt[, cnames.new], stringsAsFactors = FALSE), by = c("cell" = "cellplate"))
  # jmerge.filt <- subset(jmerge2, !is.na(louvain.facs))
  jmerge.filt <- subset(jmerge, !is.na(platename))
  
  jmerge.filt.out <- jmerge.filt
  
  # for (jcname in names(cnames.switchlist)){
  #   cname.new <- cnames.rename[[jcname]]
  #   jcname.i <- which(colnames(jmerge.filt.out) == jcname)
  #   colnames(jmerge.filt.out)[[jcname.i]] <- cname.new
  # }
  print(colnames(jmerge.filt.out))
  return(jmerge.filt.out)
  # return(jmerge.filt)
})



# Get output --------------------------------------------------------------

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/FACS_checks/H3K27me3_rep3_FACS_outputs.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

library(DescTools)

jmerge.filt.rep3 <- bind_rows(jmerge.lst) %>%
  rowwise() %>%
  mutate(CkitSca1 = Ckit ^ 2 + Sca1 ^ 2) %>%
  ungroup() %>%
  mutate(CkitSca1.win = DescTools::Winsorize(CkitSca1, probs = c(0, 0.95)),
         Lineage.win = DescTools::Winsorize(Lineage, probs = c(0.01, 0.99)))

ggplot(jmerge.filt.rep3, aes(x = CkitSca1.win)) + 
  geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.umap.proteins <- ggplot(jmerge.filt.rep3, aes(x = umap1, y = umap2, color = log10(CkitSca1.win))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K27me3, colored by Kit + Sca1") + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.umap.proteins)

m.umap.lineage <- ggplot(jmerge.filt.rep3, aes(x = umap1, y = umap2, color = log10(Lineage.win))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K27me3, colored by Lineage") + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.umap.lineage)

ggplot(jmerge.filt.rep3, aes(x = log10(Lineage), fill = cluster)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge.filt.rep3, aes(x = log10(Ckit), fill = cluster)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge.filt.rep3, aes(x = log10(Sca1), fill = cluster)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Check erythroblast  -----------------------------------------------------


jsub <- subset(jmerge.filt.rep3, cluster == "Eryths")

m0 <- ggplot(jsub, aes(x = umap1, y = umap2, color = CkitSca1.win)) + 
  geom_point() + 
  ggtitle("H3K27me3, erythroblast only") + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m1 <- ggplot(jsub, aes(x = umap1, y = umap2, color = log2(cuts_total / spikein_cuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m0.lin <- ggplot(jsub, aes(x = umap1, y = umap2, color = Lineage.win)) + 
  geom_point() + 
  ggtitle("H3K27me3, erythroblast only") + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


m2 <- ggplot(jsub, aes(x = log2(cuts_total / spikein_cuts), y = CkitSca1.win)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m2.lin <- ggplot(jsub, aes(x = log2(cuts_total / spikein_cuts), y = Lineage.win)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

multiplot(m0, m1, cols = 2)
multiplot(m0.lin, m1, cols = 2)
print(m2.lin)

print(m.umap.proteins)

ggplot(jmerge.filt.rep3, aes(x = forcats::fct_reorder(.f = cluster, .x = CkitSca1.win, .fun = median, .desc = TRUE), y = CkitSca1.win)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge.filt.rep3, aes(x = forcats::fct_reorder(.f = cluster, .x = Lineage, .fun = median, .desc = TRUE), y = Lineage)) + 
  geom_violin() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge.filt.rep3, aes(x = Ckit, y = Sca1, color = cluster)) + 
  geom_point() +
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jmerge.filt.rep3, aes(x = Lineage, fill = cluster)) + 
  geom_density(alpha = 0.25) + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~cluster) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()


# 
# 
# 
# # Test --------------------------------------------------------------------
# 
# 
# 
# 
# inf.test <- file.path(hubprefix, paste0("jyeung/data/scChiC/facs_index_data/BM_round_3/H3K4me1_H3K27me3_20200901/", jplate, "_index.csv"))
# assertthat::assert_that(file.exists(inf.test))
# 
# dat.facs <- fread(inf.test)
# 
# dat.facs <- dat.facs[!duplicated(dat.facs$Well), ]
# 
# # remove weird colnames
# cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME", "Tray X", "Tray Y", "Trigger Pulse Width",
#                    "Time 1", "Time 2", "Time 3", "Drop Phase", "ROI Bits 1-16", "ROI Bits 17-32", "Classifier Bits", "Sort Result Bits",
#                    "Sort Enable Bits")
# 
# 
# # cnames.keep <- !colnames(dat.facs) %in% cnames.remove
# cnames.keep <- cnames.sorting
# dat.facs.mat <- as.data.frame(dat.facs[, ..cnames.keep])
# cellnumbr <- sapply(dat.facs$Well, function(x) WellToCellNumber(x, totalcols = 24, shift = -1))
# rownames(dat.facs.mat) <- cellnumbr
# feats.keep <- which(apply(dat.facs.mat, 2, var) > 0)
# dat.facs.mat.filt <- scale(log2(dat.facs.mat[, feats.keep]), center = TRUE, scale = TRUE)
# 
# dat.facs.umap <- DoUmapAndLouvain(topics.mat = dat.facs.mat.filt, jsettings = jsettings)
# dat.facs.umap$cellplate <- paste(jplate, dat.facs.umap$cell, sep = "_")
# 
# ggplot(dat.facs.umap, aes(x = umap1, y = umap2)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # Merge with known cell types  --------------------------------------------
# 
# jmerge <- left_join(dat.metas$H3K27me3, dat.facs.umap %>% dplyr::rename(umap1.facs = umap1, umap2.facs = umap2, louvain.facs = louvain), by = c("cell" = "cellplate"))
# jmerge2 <- left_join(jmerge, data.frame(cellplate = paste(jplate, rownames(dat.facs.mat.filt), sep = "_"), dat.facs.mat.filt[, cnames.sorting], stringsAsFactors = FALSE), by = c("cell" = "cellplate"))
# jmerge.filt <- subset(jmerge2, !is.na(louvain.facs))
# 
# # rename?j
# jmerge.filt.out <- jmerge.filt
# 
# cnames.switchlist <- list("X..561..585.29" = "Lineage",
#                           "X..640..670.30" = "Ckit",
#                           "X..561..750.LP" = "Sca1")
# # names(cnames.switchlist) <- cnames.switchlist
# cnames.rename <- hash::hash(cnames.switchlist)
# 
# for (jcname in names(cnames.switchlist)){
#   cname.new <- cnames.rename[[jcname]]
#   jcname.i <- which(colnames(jmerge.filt.out) == jcname)
#   colnames(jmerge.filt.out)[[jcname.i]] <- cname.new
# }
# 
# return(jmerge.filt.out)
# 
# 
# 
# ggplot(jmerge.filt, aes(x = umap1.facs, y = umap2.facs, color = batch)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # ggplot(jmerge.filt, aes(x = log2(X..561..585.29 + 1), y = log2(X..640..670.30 + 1), color = batch)) + 
# #   geom_point() + 
# #   theme_bw() + 
# #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(jmerge.filt, aes(x = X..561..585.29, y = X..640..670.30, color = batch)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(jmerge.filt, aes(x = X..561..585.29, fill = batch)) + 
#   geom_density(alpha = 0.5) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(jmerge.filt, aes(x = X..561..585.29, y = X..640..670.30, color = batch)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(jmerge.filt, aes(x = X..561..750.LP, y = X..640..670.30, color = batch)) + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 


# Jake Yeung
# Date of Creation: 2021-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/FACS/explore_FACS_doubleincubate.R
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



# Load doubles ------------------------------------------------------------

inf.rdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/double_stain_outputs/scchix_outputs.H3K4me1-H3K9me3.2021-02-11.setseed.RData"
load(inf.rdata, v=T)


# Load FACS  --------------------------------------------------------------


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
jplates <- unique(subset(dat.umap.long.annot, jrep == "rep3")$experi)
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

indir.facs <- file.path(hubprefix, "jyeung/data/scChiC/facs_index_data/BM_round_3/H3K9me3_20200915")
  
jmerge.lst <- lapply(jplates, function(jplate){
  print(jplate)
  
  inf.test <- file.path(indir.facs, paste0(jplate, ".csv"))
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
  
  jmerge <- left_join(dat.umap.long.annot, dat.facs.umap, by = c("cell" = "cellplate"))
  jmerge.filt <- subset(jmerge, !is.na(platename))
  
  jmerge.filt.out <- jmerge.filt
  
  print(colnames(jmerge.filt.out))
  return(jmerge.filt.out)
})



# Get output --------------------------------------------------------------

library(DescTools)

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/FACS_checks/H3K4me1_H3K9me3_FACS_outputs.", Sys.Date(), ".pdf")

pdf(outpdf, useDingbats = FALSE)

jmerge.filt.rep3 <- bind_rows(jmerge.lst) %>%
  rowwise() %>%
  mutate(CkitSca1 = Ckit ^ 2 + Sca1 ^ 2) %>%
  ungroup() %>%
  mutate(CkitSca1.win = DescTools::Winsorize(CkitSca1, probs = c(0, 0.95)),
         Lineage.win = DescTools::Winsorize(Lineage, probs = c(0.01, 0.99)))

jmerge.filt.rep3.cleaned <- left_join(fits.merge, subset(jmerge.filt.rep3, select = c(cell, CkitSca1, Lineage, CkitSca1.win, Lineage.win, Ckit, Sca1)))

ggplot(jmerge.filt.rep3, aes(x = CkitSca1.win)) + 
  geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.umap.proteins <- ggplot(jmerge.filt.rep3, aes(x = umap1, y = umap2, color = log10(CkitSca1.win))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me1+H3K9me3, colored by Kit + Sca1") + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.umap.proteins)

# ggplot(jmerge.filt.rep3, aes(x = louv.act.impute, y = log10(CkitSca1.win))) + 
#   geom_boxplot() + 
#   geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.grid <- ggplot(jmerge.filt.rep3 %>% arrange(CkitSca1.win), aes(x = louv.act.impute, y = louv.repress.impute, color = log10(CkitSca1.win))) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

m.grid.lin <- ggplot(jmerge.filt.rep3 %>% arrange(Lineage.win), aes(x = louv.act.impute, y = louv.repress.impute, color = log10(Lineage.win))) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.lin)



m.umap.lineage <- ggplot(jmerge.filt.rep3, aes(x = umap1, y = umap2, color = log10(Lineage.win))) + 
  geom_point() + 
  theme_bw() + 
  ggtitle("H3K4me1+H3K9me3, colored by Lineage") + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.umap.lineage)

ggplot(jmerge.filt.rep3, aes(x = log10(Lineage), fill = louv.act.impute)) + 
  geom_histogram(alpha = 0.25) + 
  facet_wrap(~louv.act.impute) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()
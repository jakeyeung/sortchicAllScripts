# Jake Yeung
# Date of Creation: 2020-08-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/2-check_mouse_BM_glmpca.R
# Check mouse BM glmpca

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(glmpca)

library(hash)
library(igraph)
library(umap)

library(reshape2)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load annots BM ----------------------------------------------------------

# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
# jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jmark <- "H3K4me1"
# jdist <- "1000"
jdists <- c("0", "1000", "10000")
jsuffix <- "RemoveSmallPeaks"
names(jdists) <- paste("dist_", jdists, sep = "")


# Add activity  -----------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
# jmarks <- c("H3K4me1"); names(jmarks) <- jmarks
zscores.wide.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  mdir.lst <- lapply(jdists, function(jdist){
    print(jdist)
    mdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_output/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.", jsuffix, "-hiddenDomains_motevo_merged.closest.long.distfilt_", jdist, ".scale_0.center_0.byrow_0.", jmark, "/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.", jsuffix)
    assertthat::assert_that(dir.exists(mdir))
    return(mdir)
  })
  
  mara.out.lst <- lapply(mdir.lst, function(mdir){
    mara.out <- LoadMARA(mdir = mdir, make.cnames = FALSE)
    return(mara.out)
  })
  
  zscores.lst <- lapply(mara.out.lst, function(mara.out){
    mara.out$zscores
  })
  
  zscores.long <- lapply(names(jdists), function(jname){
    zscores <- zscores.lst[[jname]]
    zscores$jdist <- jname
    return(zscores)
  }) %>%
    bind_rows()
  zscores.long$mark <- jmark
  return(zscores.long)
  # print(head(zscores.long))
  # zscores.wide <- dcast(zscores.long, formula = motif ~ jdist, value.var = "zscore")
  # zscores.wide$mark <- jmark
  # return(zscores.wide)
})

zscores.long <- zscores.wide.lst %>%
  bind_rows() %>%
  filter(mark != "H3K27me3")

zscores.wide1 <- dcast(zscores.long %>% filter(jdist == "dist_0"), formula = motif ~ mark, value.var = "zscore") %>%
  mutate(jdist = "dist_0")
zscores.wide2 <- dcast(zscores.long %>% filter(jdist == "dist_10000"), formula = motif ~ mark, value.var = "zscore") %>%
  mutate(jdist = "dist_10000")

zscores.wide <- bind_rows(zscores.wide1, zscores.wide2)

ggplot(zscores.wide, aes(x = H3K4me1, y = H3K4me3)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~jdist)

# label relevant ones 

zscores.wide <- zscores.wide %>%
  rowwise() %>%
  mutate(motif.lab = ifelse(H3K4me1 >= 1.25 | H3K4me3 >= 1.25, motif, NA))

library(ggrepel)
ggplot(zscores.wide %>% filter(jdist == "dist_0"), aes(x = H3K4me1, y = H3K4me3, label = motif.lab)) + 
  geom_point() + 
  geom_text_repel(size = 6) + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(zscores.wide, aes(x = H3K4me1, y = H3K4me3, label = motif.lab)) + 
  geom_point() + 
  geom_text_repel(size = 6) + 
  theme_bw(12) + 
  facet_wrap(~jdist) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# # plot output
# zscores.wide.all <- zscores.wide.lst %>%
#   bind_rows()



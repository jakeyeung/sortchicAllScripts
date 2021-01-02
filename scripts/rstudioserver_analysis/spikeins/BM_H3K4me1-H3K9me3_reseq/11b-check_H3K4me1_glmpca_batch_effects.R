# Jake Yeung
# Date of Creation: 2020-12-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/11b-check_H3K4me1_glmpca_batch_effects.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# niter <- "500"
# niter2 <- "500"
# niter2 <- "1000"
# binskeep <- "0"
# binskeep2 <- "1000"
# jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
# jsuffix2 <- paste0("bincutoff_0.binskeep_", binskeep2, ".byplate.szname_none.niter_", niter2, ".reorder_rownames.dupfilt")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

dat.metas <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt"))
  fread(inf.meta)
})


binskeep <- "1000"
niter <- "1000"

# hubprefix <- "/home/jyeung/hub_oudenaarden"
# inf.k4me1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun/glmpca.H3K4me1.bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt.RData")
inf.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun.byrep/glmpca.H3K4me1.bincutoff_0.binskeep_", binskeep, ".platename_jrep.szname_none.niter_", niter, ".reorder_rownames.dupfilt.RData"))
# inf.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/same_annot_file_rerun.byrep/glmpca.H3K4me1.bincutoff_0.binskeep_", binskeep, ".platename_jrep.szname_none.niter_", niter, ".reorder_rownames.dupfilt.RData"))
print(inf.k4me1)
assertthat::assert_that(file.exists(inf.k4me1))
load(inf.k4me1, v=T)

glmfactors <- glm.out$factors

spread <- 8

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$spread <- spread
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(glmfactors, jsettings = jsettings)

dat.umap.lst <- list()
dat.umap.lst$H3K4me1 <- dat.umap
jmark <- "H3K4me1"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

m <- ggplot(dat.umap.lst[[jmark]] %>% left_join(., dat.metas[[jmark]] %>% dplyr::select(c("cell", "cluster"))), 
            aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.2, alpha = 0.2) + 
  ggtitle(jmark, paste(paste0("spread=", jsettings$spread), paste0("\nmindist=", jsettings$min_dist), sep = ",")) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# test 
mcheck <- ggplot(dat.umap.lst[[jmark]] %>% left_join(., dat.metas[[jmark]] %>% dplyr::select(c("cell", "cluster", "jrep"))), 
                 aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.8, alpha = 0.8) + 
  facet_wrap(~jrep) + 
  ggtitle(jmark, paste(paste0("spread=", jsettings$spread), paste0("\nmindist=", jsettings$min_dist), sep = ",")) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

print(m)
print(mcheck)


# Check -------------------------------------------------------------------



# Jake Yeung
# Date of Creation: 2020-02-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/9-downstream_GLMPCA_BM_KeepMorePlates.R
# description

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(glmpca)
library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmarks <- jmarks[[3]]

# jmark <- "H3K4me1"

# infs <- ""
jsuffix <- "KeepMorePlates"

# jdates <- c("2020-02-02", "2020-02-04", "2020-02-04", "2020-02-04")
# names(jdates) <- jmarks
jdate <- c("2020-02-04")

for (jmark in jmarks){
  print(jmark)
  
niter <- 1000
topn <- 150
jbins.keep <- 500
# calculating var raw
binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * mergesize


outdir <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs"
pdfdir <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.downstream"
dir.create(pdfdir)
outbase <- paste0("PZ_", jmark, ".", jsuffix, ".GLMPCA_var_correction.150.", jdate, ".binskeep_1000.devmode")
# outbase <- paste0("PZ_", jmark, ".GLMPCA_var_correction.", topn, ".", Sys.Date(), ".binskeep_", jbins.keep, ".devmode")
outpdf <- file.path(pdfdir, paste0(outbase, ".pdf"))

# outf <- "/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs/PZ_H3K4me3.KeepMorePlates.GLMPCA_var_correction.150.2020-02-04.binskeep_1000.devmode.RData"
outf <- file.path(outdir, paste0(outbase, ".RData"))

assertthat::assert_that(file.exists(outf))


load(outf, v=T)

# inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt/lda_outputs.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.AllMerged.K-30.Robj")
# inf <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs/PZ_", jmark, ".", jsuffix, ".GLMPCA_var_correction.150.2020-02-04.binskeep_1000.devmode.RData")
inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_", jmark, ".varcutoff_0.3.KeepAllPlates.K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.KeepAllPlates.K-30.Robj")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

tm.result <- posterior(out.lda)

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))

jchromos <- paste("chr", c(seq(19)), sep = "")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)


# Plot output -------------------------------------------------------------

topics.mat <- glm.out$factors

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

dat.umap.long <- dat.umap.long %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), jsep = "_")) %>%
  left_join(., dat.var)

m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) 

m.louv.plate <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)  + 
  facet_wrap(~plate)

m.var <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) 

m.var.plate <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)  + 
  facet_wrap(~plate)


pdf(outpdf,useDingbats = FALSE)
  print(m.louv)
  print(m.louv.plate)
  print(m.var)
  print(m.var.plate)
dev.off()


}

print(Sys.time() - jstart)




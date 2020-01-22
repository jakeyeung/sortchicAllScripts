# Jake Yeung
# Date of Creation: 2020-01-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/6-correct_variance_LDA_k4me3.R
# K4me3



library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(JFuncs)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)



# Constants ---------------------------------------------------------------

binsize <- 50000
mergesize <- 1000
bigbinsize <- 50000 * 1000

jsystem <- "BoneMarrow"
jmark <- "H3K4me3"

jcutoff.ncuts.var <- 0.3
# bad.plates.grep <- paste("Bl6-BM-Linneg-H3K4me1", "ChIC-Bl6-BM-H3K4me1", sep = "|")
bad.plates.grep <- ""

outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/quality_control_postLDA_var_BM"
outname <- paste0("BM_", jmark, ".varcutoff_", jcutoff.ncuts.var, ".KeepAllPlates.rds")
pdfname <- paste0("BM_", jmark, ".varcutoff_", jcutoff.ncuts.var, ".KeepAllPlates.pdf")
outf <- file.path(outdir, outname)
outpdf <- file.path(outdir, pdfname)

# assertthat::assert_that(!file.exists(outf))
# assertthat::assert_that(!file.exists(outpdf))


# Set up  -----------------------------------------------------------------


inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-12.bsizestepsize_50000_50000.NoSliding/lda_outputs.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoff_1000-500-1000-1000.TAcutoff_0.5.K-30.Robj")
load(inf, v=T)


# Show UMAP ---------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


tm.result <- posterior(out.lda)

# topic 7 contains intrachrom var? 
# tm.result$topics <- tm.result$topics[, -7]
# tm.result$terms <- tm.result$terms[-7, ]

topics.mat <- tm.result$topics


umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.umap <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)


dat.merge <- left_join(dat.umap.long, dat.var) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "-"),
         plate = ClipLast(cell, jsep = "_"),
         # prefix = gsub("PZ-Bl6-BM", "Linneg", paste(strsplit(gsub("PZ-", "", cell), "-")[[1]][1:4], collapse = "-")))
         prefix = paste(strsplit(cell, "-")[[1]][1:4], collapse = "-"))

print(unique(dat.merge$prefix))
print(unique(dat.merge$experi))

m.umap.plates <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~plate)

m.umap.prefix <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) + 
  facet_wrap(~prefix)

m.umap.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)

m.umap.var.plates <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~plate)

# Show raw  ---------------------------------------------------------------


dat.var.raw <- CalculateVarRaw(as.matrix(count.mat), merge.size = mergesize, chromo.exclude.grep = "^chrX|^chrY", jpseudocount = 1, jscale = 10^6, calculate.ncuts = TRUE)

dat.merge2 <- left_join(dat.merge, dat.var.raw)

# Correlate raw intrachrom var with UMAP  ---------------------------------

m.rawvar.vs.imputevar <- ggplot(dat.merge2, aes(x = ncuts.var, y = cell.var.within.sum.norm, color = prefix)) + geom_point(alpha = 0.5) + 
  scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste("Intrachromo var raw from", bigbinsize/10^6, "MB bins")) +  
  ylab("Imputed intrachromo var from LDA") + 
  ggtitle(jmark, jsystem)
print(m.rawvar.vs.imputevar)

m.rawvar.vs.imputevar.plates <- ggplot(dat.merge2, aes(x = ncuts.var, y = cell.var.within.sum.norm, color = prefix)) + geom_point(alpha = 0.5) + 
  scale_x_log10() + scale_y_log10() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(paste("Intrachromo var raw from", bigbinsize/10^6, "MB bins")) +  
  ylab("Imputed intrachromo var from LDA") + facet_wrap(~plate) + 
  ggtitle(jmark, jsystem)
print(m.rawvar.vs.imputevar.plates)

# show ncuts vs var
m.ncutsVSvar <- ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var, color = prefix)) + geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10()  

m.ncutsVSvar.plate <- ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var, color = prefix)) + geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + facet_wrap(~plate)




# filter bad plates and bad cells and redo? 
unique(dat.merge2$prefix)

# calculate intrachromo cutoff 

dat.linneg <- subset(dat.merge2, grepl("Linneg", plate))

m.var.cutoff <- NULL

m.ncutsVSvar.plate.cutoff <- ggplot(dat.merge2, aes(x = ncuts, y = ncuts.var, color = prefix)) + geom_point(alpha = 0.25)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() + facet_wrap(~plate) + geom_hline(yintercept = jcutoff.ncuts.var)


# Write new count mat -----------------------------------------------------


if (bad.plates.grep != ""){
  dat.keep <- subset(dat.merge2, !grepl(bad.plates.grep, plate) & ncuts.var > jcutoff.ncuts.var)
} else {
  dat.keep <- subset(dat.merge2, ncuts.var > jcutoff.ncuts.var)
}
cells.keep <- dat.keep$cell

print(unique(dat.keep$prefix))

count.mat.keep <- count.mat[, cells.keep]

print(dim(count.mat))
print(dim(count.mat.keep))

saveRDS(count.mat.keep, outf)





# Correct LDA -------------------------------------------------------------

# topics.mat.filt <- topics.mat[cells.keep, ]
topics.mat.input <- scale(topics.mat, center = TRUE, scale = TRUE)

pca.out <- prcomp(topics.mat.input, center = FALSE, scale. = FALSE)

print(summary(pca.out))

dat.pca <- data.frame(cell = rownames(pca.out$x), pca.out$x, stringsAsFactors = FALSE) %>%
  left_join(., dat.merge2)

m.pca <- ggplot(dat.pca, mapping = aes(x = PC1, y = PC2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
m.pca.var <- ggplot(dat.pca, mapping = aes(x = PC1, y = PC2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



rot.mat <- pca.out$rotation[, 2:15]

topics.mat.rot <- topics.mat.input %*% rot.mat

umap.out <- umap(topics.mat.rot, config = jsettings)
dat.umap.long.rot <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long.rot <- DoLouvain(topics.mat.rot, custom.settings.louv = jsettings, dat.umap.long.rot, clstr.cname = "new_louvain")
dat.umap.long.rot <- dat.umap.long.rot %>%
  left_join(., subset(dat.merge2, select = c(-umap1, -umap2)), by = c("cell"))


m.umap.cor <- ggplot(dat.umap.long.rot, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

m.umap.cor.var <- ggplot(dat.umap.long.rot, aes(x = umap1, y = umap2, color = cell.var.within.sum)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_viridis_c(direction = -1)

m.umap.cor.plates <- ggplot(dat.umap.long.rot, aes(x = umap1, y = umap2, color = prefix)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) +
  facet_wrap(~plate)

m.umap.cor.prefix <- ggplot(dat.umap.long.rot, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) +
  facet_wrap(~prefix)

multiplot(m.umap.prefix, m.umap.cor.prefix, cols = 1)

m.umap.cor.newlouvain <- ggplot(dat.umap.long.rot, aes(x = umap1, y = umap2, color = new_louvain)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette) +
  facet_wrap(~prefix)






# Write plots -------------------------------------------------------------

pdf(file = outpdf, useDingbats = FALSE)

  m.umap
  m.umap.var
  m.umap.var.plates

  m.ncutsVSvar
  m.ncutsVSvar.plate
  m.ncutsVSvar.plate.cutoff

  m.rawvar.vs.imputevar
  m.rawvar.vs.imputevar.plates
  m.var.cutoff
  
  plot(pca.out$sdev ^ 2 / sum(pca.out$sdev ^ 2))
  
  m.pca
  m.pca.var

  m.umap.cor
  m.umap.cor.var
  m.umap.cor.plates
  m.umap.cor.prefix
  m.umap.cor.newlouvain

  multiplot(m.umap.prefix, m.umap.cor.prefix, cols = 1)


dev.off()







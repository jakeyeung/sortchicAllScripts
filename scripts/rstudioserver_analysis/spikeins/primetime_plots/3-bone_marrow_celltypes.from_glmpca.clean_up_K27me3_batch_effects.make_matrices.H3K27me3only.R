# Jake Yeung
# Date of Creation: 2020-12-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_celltypes.from_glmpca.clean_up_K27me3_batch_effects.make_matrices.H3K27me3only.R
# description
rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)
library(scchicFuncs)
library(JFuncs)

jstart <- Sys.Date()


# Load objs  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmark.test <- "H3K4me3"
# jmark.test <- "H3K4me1"
jmark.test <- "H3K27me3"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2"
dir.create(outdir)
outrdata <- file.path(outdir, paste0("matrix_imputed_mark_", jmark.test, ".RData"))
outpdf <- file.path(outdir, paste0("matrix_imputed_mark_", jmark.test, ".pdf"))

inf.lda <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep1rep2rep3reseq_varfilt_2020-12-15/lda_outputs.mat_H3K27me3_rep1rep2rep3reseq.TSS.K-30.binarize.FALSE/ldaOut.mat_H3K27me3_rep1rep2rep3reseq.TSS.K-30.Robj"
load(inf.lda, v=T)

# get imputed
tm.result <- posterior(out.lda)
dat.imputed.log2 <- log2(t(tm.result$topics %*% tm.result$terms))

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.metadata.txt"
dat.meta <- fread(inf.meta)

pdf(outpdf, useDingbats = FALSE)


# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

tm.result.lst <- list()
tm.result.lst[[jmark.test]] <- tm.result

count.mat.raw.lst <- list()
count.mat.raw.lst[[jmark.test]] <- count.mat

inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
dat.genes <- fread(inf.genes)

dat.meta.lst <- list()
dat.meta.lst[[jmark.test]] <- dat.meta

# dat.meta.lst <- lapply(jmarks, function(jmark){
#   inf.meta <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/metadata_umap_celltype_cuts.", jmark, ".txt")
#   fread(inf.meta)
# })

dat.imputed.log.lst <- lapply(tm.result.lst, function(tm.result){
  log2(t(tm.result$topics %*% tm.result$terms))
})

dat.imputed.log.lst <- list()
dat.imputed.log.lst[[jmark.test]] <- dat.imputed.log2


# Make heatmap, fix batch effects -----------------------------------------

cells.keep <- dat.meta.lst[[jmark.test]]$cell

genes.keep <- dat.genes$gene

# rkeep <- genes.keep %in% rownames()

mat4hm <- dat.imputed.log.lst[[jmark.test]][genes.keep, cells.keep]
mat4hm.uniq <- dat.imputed.log.lst[[jmark.test]][unique(genes.keep), cells.keep]

heatmap3::heatmap3(mat4hm, Rowv = NA, Colv = NA, ColSideColors = dat.meta.lst[[jmark.test]]$colorcode, scale = "row", RowSideColors = dat.genes$colorcode, revC = TRUE, main = paste0(jmark.test))


# Fix batch effects?  -----------------------------------------------------


mat4hm.long <- data.table::melt(mat4hm.uniq)
colnames(mat4hm.long) <- c("rname", "cell", "log2exprs")

# mat4hm.long <- subset(mat4hm.long, rname %in% unique(genes.keep))

mat4hm.long.annot <- left_join(mat4hm.long, dat.meta.lst[[jmark.test]])

mat4hm.long.annot$jrep2 <- sapply(mat4hm.long.annot$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))

jgene <- "Hbb"
jgene <- "Hbb"

jgene <- "Hlf"

jgene <- "Tal1"
jgene <- "S100a8"
jgene <- "Cebpb"
jgene <- "Gcnt2"

jgene <- "Pax5"

jgene <- "S100a8"
jrname <- grep(pattern = jgene, x = mat4hm.long.annot$rname, value = TRUE)[[1]]
print(jrname)
jfilt <- subset(mat4hm.long.annot, rname == jrname)

ggplot(jfilt, aes(x = cluster, y = log2exprs)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~jrep) + 
  ggtitle(jmark.test, jrname) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Fix the batch effect  ---------------------------------------------------


# jfilt <- jfilt %>%
#   rowwise() %>%
#   mutate(jrep2 = ifelse(jrep == "rep1old", "zold", "anew"),
#          jrep3 = jrep,
#          jrep3 = gsub("rep3", "arep3old", jrep3),
#          jrep3 = gsub("rep2", "brep2old", jrep3),
#          jrep3 = gsub("rep1", "crep1old", jrep3))
  
jout <- lm(data = jfilt, formula = log2exprs ~ 1 + jrep2 + cluster + jrep2:cluster)
# jout <- lm(data = jfilt, formula = log2exprs ~ 1 + jrep3 + cluster + jrep3:cluster)

# effects to back out 
jout$coefficients[["jrep2zold:clusterBcells"]]

plateeffect.i <- grep("jrep2.*.old$", names(jout$coefficients), value = FALSE)
plateeffect <- jout$coefficients[plateeffect.i]
interactioneffect.i <- grep("jrep2zold:", names(jout$coefficients), value = FALSE)
interactioneffect <- jout$coefficients[interactioneffect.i]


jfilt <- jfilt %>%
  rowwise() %>%
  mutate(plateadj2 = AssignEffect(jrep2, cluster, plateeffect, jtype = "plate"),
         clstradj2 = AssignEffect(jrep2, cluster, interactioneffect, jtype = "cluster"))

jfilt <- jfilt %>%
  rowwise() %>%
  mutate(plateadj = ifelse(jrep == "rep1old", plateeffect, 0),
         clstradj = ifelse(jrep == "rep1old" & cluster != "Basophils", interactioneffect[[paste0("jrep2zold:cluster", cluster)]], 0))

ggplot(jfilt, aes(x = cluster, y = log2exprs - plateadj - clstradj)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~jrep) + 
  ggtitle(jmark.test, jrname) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jfilt, aes(x = cluster, y = log2exprs - plateadj2 - clstradj2)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~jrep) + 
  ggtitle(jmark.test, jrname) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfilt.adj <- AdjustBatchEffect(jfilt, plateprefix = "jrep2", clusterprefix = "cluster", platesuffix = "old")

ggplot(jfilt.adj, aes(x = cluster, y = log2exprs - plateadj - clstradj)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~jrep) + 
  ggtitle(jmark.test, jrname) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


ggplot(jfilt.adj, aes(x = cluster, y = log2exprs - plateadj - clstradj)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  facet_wrap(~jrep) + 
  ggtitle(jmark.test, jrname) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# Do for all  -------------------------------------------------------------

dat.adj <- mat4hm.long.annot %>%
  group_by(rname) %>%
  do(AdjustBatchEffect(.)) 

dat.adj <- dat.adj %>%
  rowwise() %>%
  mutate(log2exprsadj = log2exprs - plateadj2 - clstradj2)

# rnames.unique <- as.character(unique(dat.adj$rname))
# dat.adj.filt <- subset(dat.adj, rname %in% rnames.unique)

mat.adj <- data.table::dcast(dat.adj, formula = rname ~ cell, value.var = "log2exprsadj")

save(dat.adj, mat.adj, file = outrdata)

# remake mat
mat4hm.adj <- as.data.frame(mat.adj)
rownames(mat4hm.adj) <- mat4hm.adj$rname
mat4hm.adj$rname <- NULL

mat4hm.adj <- mat.adj[genes.keep, cells.keep]
heatmap3::heatmap3(mat4hm.adj, Rowv = NA, Colv = NA, ColSideColors = dat.meta.lst[[jmark.test]]$colorcode, scale = "row", RowSideColors = dat.genes$colorcode, revC = TRUE, main = paste0(jmark.test, "adjusted"))

dev.off()

print(Sys.Date() - jstart)

# Plot heatmap adjusted ---------------------------------------------------





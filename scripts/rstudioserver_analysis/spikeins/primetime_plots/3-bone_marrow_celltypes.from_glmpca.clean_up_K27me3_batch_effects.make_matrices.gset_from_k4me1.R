# Jake Yeung
# Date of Creation: 2020-12-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_celltypes.from_glmpca.clean_up_K27me3_batch_effects.make_matrices.gset_from_k4me1.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

jstart <- Sys.Date()


AssignEffect <- function(batch, cluster, effects, jtype = c("plate", "cluster"), plateprefix = "jrep2", clusterprefix = "cluster"){
  if (jtype == "plate"){
    jname <- paste0(plateprefix, batch)
    if (jname %in% names(effects)){
      # print("Adjusting plate nonzero")
      indx <- names(effects) %in% jname
      adj <- effects[indx]
      # adj <-  effects[[effects[jname]]]
    }  else {
      # print("Adjusting plate  zero")
      adj <- 0
    }
  } else if (jtype == "cluster"){
    jname <- paste0(plateprefix, batch, ":", clusterprefix,  cluster)
    if (jname %in% names(effects)){
      # print(jname)
      # print(names(effects))
      indx <- names(effects) %in% jname
      adj <- effects[indx]
      # print("Adjusting cluster nonzero")
    } else {
      # print("Adjusting cluster zero")
      adj <- 0
    }
  } else {
    print(paste("jtype must be plate or cluster", jtype))
  }
  return(adj)
}

AdjustBatchEffect <- function(jdat, plateprefix = "jrep2", clusterprefix = "cluster", platesuffix = "old"){
  jout <- lm(data = jdat, formula = log2exprs ~ 1 + jrep2 + cluster + jrep2:cluster)
  plateeffect.i <- grep(paste0(plateprefix, ".*.", platesuffix, "$"), names(jout$coefficients), value = FALSE)
  plateeffect <- jout$coefficients[plateeffect.i]
  # interactioneffect.i <- grep(plateprefix, ".*.:", names(jout$coefficients), value = FALSE)
  interactioneffect.i <- grep(paste0(plateprefix, ".*.:", clusterprefix), names(jout$coefficients), value = FALSE)
  interactioneffect <- jout$coefficients[interactioneffect.i]
  
  jdat <- jdat %>%
    rowwise() %>%
    mutate(plateadj2 = AssignEffect(jrep2, cluster, plateeffect, jtype = "plate"),
           clstradj2 = AssignEffect(jrep2, cluster, interactioneffect, jtype = "cluster"),
           log2exprsadj = log2exprs - plateadj2 - clstradj2)
  return(jdat)
}



# Load objs  --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark.test <- "H3K4me3"
# jmark.test <- "H3K4me1"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2.H3K4me1_TSS_topics2"
dir.create(outdir)
outrdata <- file.path(outdir, paste0("matrix_imputed_mark_", jmark.test, ".RData"))
outpdf <- file.path(outdir, paste0("matrix_imputed_mark_", jmark.test, ".pdf"))

pdf(outpdf, useDingbats = FALSE)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

inf.tmresult <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/imputed_TSS_three_marks.rds")
tm.result.lst <- readRDS(inf.tmresult)

inf.raw <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/raw_cuts_TSS_three_marks.rds")
count.mat.raw.lst <- inf.raw

# inf.genes.old <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
# inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics.txt")
# dat.genes.old <- fread(inf.genes.old)
# dat.genes <- fread(inf.genes)
# 
# # filter out and add colors
# jset2color <- hash::hash(unique(dat.genes.old$jset), unique(dat.genes.old$colorcode))
# dat.genes$colorcode <- sapply(dat.genes$jset, function(x) AssignHash(x = x, jhash = jset2color, null.fill = NA))
# dat.genes.sub <- subset(dat.genes, !is.na(colorcode))
# outf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics.filt.colorcoded.txt")
# fwrite(x = dat.genes.sub, outf.genes)

inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt")
dat.genes <- fread(inf.genes)


dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/metadata_umap_celltype_cuts.", jmark, ".txt")
  fread(inf.meta)
})

dat.imputed.log.lst <- lapply(tm.result.lst, function(tm.result){
  log2(t(tm.result$topics %*% tm.result$terms))
})


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

jgene <- "Gcnt2"
jgene <- "Hbb"
jgene <- "Hbb"
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

mat4hm.adj <- mat4hm.adj[genes.keep, cells.keep]
heatmap3::heatmap3(mat4hm.adj, Rowv = NA, Colv = NA, ColSideColors = dat.meta.lst[[jmark.test]]$colorcode, scale = "row", RowSideColors = dat.genes$colorcode, revC = TRUE, main = paste0(jmark.test, "adjusted"))

dev.off()

print(Sys.Date() - jstart)

# Plot heatmap adjusted ---------------------------------------------------





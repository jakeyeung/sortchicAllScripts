# Jake Yeung
# Date of Creation: 2021-01-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/18-analyze_scchix_model_output_tweaks.mereg_two_fits.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(ggforce)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

ctypes.order.k4 <- c("Eryths", "Bcells", "NKs", "DCs", "Granulocytes", "Basophils", "pDCs", "HSPCs")
names(ctypes.order.k4) <- ctypes.order.k4
ctypes.order.k9 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")
names(ctypes.order.k9) <- ctypes.order.k9

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/double_stain_outputs"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("scchix_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".setseed.pretty.pdf"))
outf <- file.path(outdir, paste0("scchix_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".setseed.txt"))
outrdata <- file.path(outdir, paste0("scchix_outputs.H3K4me1-H3K9me3.", Sys.Date(), ".setseed.RData"))

set.seed(0)
make.plots <- TRUE
write.tables <- TRUE

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

# Load metas  -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})

batch2col <- hash::hash(dat.metas$H3K4me1$batch, dat.metas$H3K4me1$stypecol)

# Load output -------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

wmin <- 0.1
wmax <- 0.9
pcount <- 0

wmax.vec <- c(0.5, 0.9); names(wmax.vec) <- wmax.vec

# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.4.wmax_0.6.RData")
# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.wmax_0.7.RData")
# infrdata <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_0.wmax_0.7.pseudocount_0.RData")

infsrdata <- lapply(wmax.vec, function(wmax){
  print(wmax)
  infrdata <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/scchix_outputs.constrain_weights/unmix_scchix_inputs_clstr_by_celltype_H3K4me1xH3K9me3.removeNA_FALSE.wmin_", wmin, ".wmax_", wmax, ".pseudocount_", pcount, ".RData"))
  assertthat::assert_that(file.exists(infrdata))
  return(infrdata)
})

fits.out.lst <- lapply(infsrdata, function(infrdata){
  load(infrdata, v=T)
  fits.out <- act.repress.coord.lst
  w.lst <- sapply(fits.out, function(x) x$w)
  
  # if louvains are now from clusters need eto rethink jcoord
  cell.vec <- names(fits.out)
  names(cell.vec) <- cell.vec
  coords.dbl <- lapply(cell.vec, function(jcell){
    jfit <- fits.out[[jcell]]
    jweight <- fits.out[[jcell]]$w
    p.mat <- SoftMax(jfit$ll.mat)
    jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
    jmax <- max(p.mat)
    
    # rows are active, columns are repress I THINK?
    # TODO: assumes underscores be careful!
    jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
    jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
    
    if (grepl("_", jlouv.act)){
      jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
    }
    if (grepl("_", jlouv.repress)){
      jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
    }
    out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
    return(out.dat)
  }) %>%
    bind_rows()
  return(coords.dbl)
})



# Check erythroblast estimates for both models  ---------------------------

jctype <- "Basophils"
jctype <- "pDCs"
jctype <- "HSPCs"
jctype <- "Eryths"

jctypes <- unique(fits.out.lst[[1]]$louv.act)
for (jctype in jctypes){
  dat.ctype1 <- subset(fits.out.lst[[1]], louv.act == jctype)
  dat.ctype2 <- subset(fits.out.lst[[2]], louv.act == jctype)
  print(jctype)
  print(table(dat.ctype1$louv.repress))
  print(table(dat.ctype2$louv.repress))
}


# Mix in fits  ------------------------------------------------------------


fits.eryth <- subset(fits.out.lst[["0.5"]], louv.act == "Eryths") %>%
  mutate(wmax = 0.5)
fits.noneryth <- subset(fits.out.lst[["0.9"]], louv.act != "Eryths") %>%
  mutate(wmax = 0.9)

fits.merge <- rbind(fits.noneryth, fits.eryth)

fits.merge$louv.act <- factor(fits.merge$louv.act, levels = ctypes.order.k4)
fits.merge$louv.repress <- factor(fits.merge$louv.repress, levels = ctypes.order.k9)


m.grid <- ggplot(fits.merge, aes(x = louv.act, y = louv.repress, louv.act)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)

m.grid <- ggplot(fits.merge %>% filter(lnprob == 0), aes(x = louv.act, y = louv.repress, louv.act, color = lnprob)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)



# Show UMAP  --------------------------------------------------------------

# load UMAP 
# jdate <- "2021-01-31"
# inf.umap <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/lda_and_clusters/ClusterAnnot.lda_and_datmerged.k4_k9_", jsuffix, ".H3K4me1xH3K9me3.", jdate, ".RData")
# assertthat::assert_that(file.exists(inf.umap))

inf.umap <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/lda_and_clusters/ClusterAnnot.lda_and_datmerged.k4_k9_dynamic_regions.H3K4me1xH3K9me3.2021-01-30.RData"
load(inf.umap, v=T)

dat.merge.annot <- dat.merge %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"),
         plate = 13,  # set 13 so it assigns LineageNeg or Unenriched
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         jrep = GetRepBM(experiname = experi),
         batch = AnnotateSortFromLayoutBMall(plate, rowcoord, colcoord, jrep, jmark)) %>%
  mutate(batch = ifelse(batch == "LinNeg", "Linneg", batch),
         stypecol = AssignHash(batch, batch2col, null.fill = batch))

ggplot(dat.merge.annot, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.dbl <- dat.merge.annot

dat.umap.dbl.merge <- left_join(dat.umap.dbl, fits.merge) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.check.repress <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

m.check.act <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = louv.act)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

m.check.w <- ggplot(dat.umap.dbl.merge, aes(x = umap1, y  = umap2, color = w)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_viridis_c()

JFuncs::multiplot(m.check.repress, m.check.act, cols = 2)


# Check total cuts --------------------------------------------------------

cuts.total <- data.frame(cell = colnames(count.mat), cuts = colSums(count.mat), stringsAsFactors = FALSE)

dat.umap.dbl.merge2 <- left_join(dat.umap.dbl.merge, cuts.total)

ggplot(dat.umap.dbl.merge2, aes(x = umap1, y = umap2, color = log2(cuts))) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.dbl.merge2, aes(x = umap1, y = umap2, color = w)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(density(log2(dat.umap.dbl.merge2$cuts)))

jcutoff <- 2^8
lncutoff <- -Inf
# lncutoff <- -0.01
wcutoff <- 0.7

m.check.repress2  <- ggplot(dat.umap.dbl.merge2 %>% mutate(louv.repress = ifelse(cuts > jcutoff & lnprob > lncutoff, as.character(louv.repress), NA)), 
                            aes(x = umap1, y  = umap2, color = louv.repress)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

m.check.act2 <- ggplot(dat.umap.dbl.merge2 %>% mutate(louv.repress = ifelse(cuts > jcutoff & lnprob > lncutoff, as.character(louv.repress), NA)), aes(x = umap1, y  = umap2, color = louv.act)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette, na.value = "grey85") 

JFuncs::multiplot(m.check.repress2, m.check.act2, cols = 2)



ggplot(dat.umap.dbl.merge2 %>% filter(cuts > jcutoff & lnprob > lncutoff & w < wcutoff), aes(x = umap1, y  = umap2, color = louv.repress)) + 
  geom_point() + theme_bw() + 
  scale_color_manual(values = cbPalette, na.value = "grey85")  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.grid2 <- ggplot(dat.umap.dbl.merge2 %>% filter(cuts > jcutoff & lnprob > lncutoff), aes(x = louv.act, y = louv.repress, louv.act, color = lnprob)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")

print(m.grid2)


m.grid2 <- ggplot(dat.umap.dbl.merge2, aes(x = louv.act, y = louv.repress, color = log2(cuts))) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid2)



# Get nearest neighbors ---------------------------------------------------

tm.result <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8

topics.mat <- tm.result$topics
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)

# get the nearest neighbors out

dat.umap.long.annot <- left_join(dat.umap.long, subset(dat.umap.dbl.merge2, select = c(-umap1, -umap2, -louvain)))

m.batch <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = stypecol)) + 
  geom_point(size = 2.5, alpha = 0.5) + 
  scale_color_identity() + 
  theme_minimal() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.batch)
multiplot(m.batch, m.batch, cols = 2)

m1.act <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m1.repress <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louv.repress)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m2 <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

multiplot(m1.act, m2, cols = 2)
multiplot(m1.repress, m2, cols = 2)
multiplot(m1.act, m1.repress, m2, cols = 3)



# Merge louvains  ---------------------------------------------------------


print(m2)

louvs.merge <- list("3" = "3,6",
                    "6" = "3,6",
                    "2" = "2,8",
                    "8" = "2,8")

dat.umap.long.annot$louvainmerged <- sapply(dat.umap.long.annot$louvain, function(x) ifelse(!is.null(louvs.merge[[x]]), louvs.merge[[x]], x))

m2.after <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louvainmerged)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

print(m2.after)

# fill in louv.repress
dat.umap.long.annot.repress.merge <- dat.umap.long.annot %>%
  group_by(louv.act, louv.repress) %>%
  summarise(ncell = length(louv.repress)) %>%
  group_by(louv.act) %>%
  mutate(nfrac = ncell / sum(ncell))

cell2repress <- hash::hash(dat.umap.long.annot$cell, dat.umap.long.annot$louv.repress)
cell2act <- hash::hash(dat.umap.long.annot$cell, dat.umap.long.annot$louv.act)

knn.out.repress <- umap.out$knn$indexes
rnames.orig.repress <- rownames(knn.out.repress); names(rnames.orig.repress) <- rnames.orig.repress
rownames(knn.out.repress) <- sapply(rownames(knn.out.repress), function(x) AssignHash(x, cell2repress, null.fill = x))

knn.out.act <- umap.out$knn$indexes
rnames.orig.act <- rownames(knn.out.act); names(rnames.orig.act) <- rnames.orig.act
rownames(knn.out.act) <- sapply(rownames(knn.out.act), function(x) AssignHash(x, cell2act, null.fill = x))


# check an eryth 
jsub <- subset(dat.umap.long.annot, louv.act == "Eryths" & louv.repress != "Eryths")
jsub <- subset(dat.umap.long.annot, louv.act == "NKs" & louv.repress != "Bcells")
# jsub <- subset(dat.umap.long.annot, louv.act == "NKs" & louv.repress == "Bcells")
jcell <- jsub$cell[[2]]

indx <- which(rnames.orig.repress == jcell)

sort(table(rownames(knn.out.repress)[knn.out.repress[indx, ]]), decreasing = TRUE)

# Update repress assignments ----------------------------------------------

ctype.neighbors.counts.repress <- sapply(rnames.orig.repress, function(jcell){
  indx <- which(rnames.orig.repress == jcell)
  xsummary <- sort(table(rownames(knn.out.repress)[knn.out.repress[indx, ]]), decreasing = TRUE)
  return(names(xsummary)[[1]])
})

ctype.neighbors.counts.act <- sapply(rnames.orig.act, function(jcell){
  indx <- which(rnames.orig.act == jcell)
  xsummary <- sort(table(rownames(knn.out.act)[knn.out.act[indx, ]]), decreasing = TRUE)
  return(names(xsummary)[[1]])
})

dat.umap.long.annot$louv.repress.impute <- sapply(dat.umap.long.annot$cell, function(x) ctype.neighbors.counts.repress[[x]])
dat.umap.long.annot$louv.act.impute <- sapply(dat.umap.long.annot$cell, function(x) ctype.neighbors.counts.act[[x]])

ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louv.repress.impute)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louv.act.impute)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = louv.act)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


m.grid.impute <- ggplot(dat.umap.long.annot, aes(x = louv.act.impute, y = louv.repress.impute, color = log2(cuts))) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_viridis_c() + 
  theme(aspect.ratio=0.6) +
  xlab("Active Clusters") + 
  ylab("Repressive Clusters") + 
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.impute)

# color by lineage

m.grid.impute.cellsurface <- ggplot(dat.umap.long.annot, aes(x = louv.act.impute, y = louv.repress.impute, color = stypecol)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  # scale_color_viridis_c() + 
  scale_color_identity() + 
  theme(aspect.ratio=0.6) +
  xlab("Active Clusters") + 
  ylab("Repressive Clusters") + 
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.impute.cellsurface)

# calculate fractions
dat.umap.long.annot.frac <- dat.umap.long.annot %>%
  group_by(louv.act.impute, stypecol, batch) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louv.act.impute) %>%
  mutate(fraccells = ncells / sum(ncells))

dat.umap.long.annot.frac.repress <- dat.umap.long.annot %>%
  group_by(louv.repress.impute, stypecol, batch) %>%
  summarise(ncells = length(cell)) %>%
  group_by(louv.repress.impute) %>%
  mutate(fraccells = ncells / sum(ncells))

ggplot(dat.umap.long.annot.frac, aes(x = louv.act.impute, y = fraccells, fill = stypecol)) + 
  geom_col() + 
  ggtitle("louv active") + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.umap.long.annot.frac.repress, aes(x = louv.repress.impute, y = fraccells, fill = stypecol)) + 
  geom_col() + 
  ggtitle("louv repress") + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# color by celltype 

cluster2col.repress <- hash::hash(dat.metas$H3K9me3$cluster, dat.metas$H3K9me3$clustercol)
cluster2col.act <- hash::hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
dat.umap.long.annot$clustercol.repress <- sapply(as.character(dat.umap.long.annot$louv.repress.impute), function(x) cluster2col.repress[[x]])
dat.umap.long.annot$clustercol.act <- sapply(as.character(dat.umap.long.annot$louv.act.impute), function(x) cluster2col.act[[x]])

# arrange clusters

ctypes.act <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "DCs", "pDCs", "HSPCs"); names(ctypes.act) <- ctypes.act
ctypes.repress <- c("Eryths", "Bcells", "Granulocytes", "HSPCs"); names(ctypes.repress) <- ctypes.repress

dat.umap.long.annot$louv.act.impute <- factor(as.character(dat.umap.long.annot$louv.act.impute), levels = ctypes.act)
dat.umap.long.annot$louv.repress.impute <- factor(as.character(dat.umap.long.annot$louv.repress.impute), levels = ctypes.repress)

m.grid.impute <- ggplot(dat.umap.long.annot, aes(x = louv.act.impute, y = louv.repress.impute, color = clustercol.repress)) +
  geom_point(alpha = 0.25, position = position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=0.6, legend.position = "bottom") +
  xlab("H3K4me1 Clusters") + 
  ylab("H3K9me3 Clusters") + 
  scale_color_identity(guide = "legend",
                       labels = unique(dat.umap.long.annot$louv.repress.impute), 
                       breaks = unique(dat.umap.long.annot$clustercol.repress)) +
  # scale_color_manual(values = cbPalette) + xlab(paste0(jmarks[[1]], " Clusters (arbitrary names)")) + ylab(paste0(jmarks[[2]], " Clusters (arbitrary names)")) +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.impute)

# show UMAPs 
m1.act <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = clustercol.act)) + 
  geom_point(size = 2.5) + 
  theme_minimal() + 
  scale_color_identity(guide = "legend",
                       labels = unique(dat.umap.long.annot$louv.act.impute),
                       breaks = unique(dat.umap.long.annot$clustercol.act)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m1.repress <- ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = clustercol.repress)) + 
  geom_point(size = 2.5) + 
  theme_minimal() + 
  scale_color_identity(guide = "legend",
                       labels = unique(dat.umap.long.annot$louv.repress.impute), 
                       breaks = unique(dat.umap.long.annot$clustercol.repress)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

print(m1.act)
print(m1.repress)

multiplot(m1.act, m1.repress, cols = 2)


# Plot highlighting different cells  --------------------------------------

jcells <- c("PZ-BM-rep3-H3K9me3-H3K4me1-4_179", "PZ-BM-rep3-H3K9me3-H3K4me1-1_290", "PZ-BM-rep3-H3K9me3-H3K4me1-1_133", "PZ-BM-rep3-H3K9me3-H3K4me1-1_46", "PZ-BM-rep3-H3K9me3-H3K4me1-1_140", "PZ-BM-rep3-H3K9me3-H3K4me1-2_236", "PZ-BM-rep3-H3K9me3-H3K4me1-4_352", "PZ-BM-rep3-H3K9me3-H3K4me1-1_8")
names(jcells) <- jcells

for (jcell in jcells){
  mcell <- ggplot(dat.umap.long.annot %>% rowwise() %>% mutate(is.cell = cell == jcell) %>% arrange(is.cell), aes(x = umap1, y = umap2, color = is.cell, size = is.cell)) + 
    geom_point() + 
    ggtitle(jcell) + 
    theme_minimal() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(mcell)
}

if (make.plots){
  dev.off()
}


# Write tables output -----------------------------------------------------

if (write.tables){
  fwrite(dat.umap.long.annot, file = outf, sep = "\t")
  save(fits.merge, dat.umap.long.annot, out.lda, file = outrdata)
}

# Jake Yeung
# Date of Creation: 2019-06-07
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/make_scChIC_signal_across_chromosome.R
# Show scChIC signal across chromosome 



library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scran)
# try CONOS
library(conos)
library(Seurat)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Constants ---------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"
outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/trajectories_stringent"

# Load LDA objects --------------------------------------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
infs.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"

infs$H3K4me3 <- infs.stringent

tm.result.lst <- lapply(infs, LoadGetTmResult)



# Load trajectories -------------------------------------------------------

inf.traj.stringent <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.traj.stringent, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs$H3K4me3
trajs.stringent <- trajs$H3K4me3
trajs.objs.stringent <- trajs.objs$H3K4me3

inf.traj <- "/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata"
load(inf.traj, v=T)

dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent
trajs$H3K4me3 <- trajs.stringent
trajs.objs$H3K4me3 <- trajs.objs.stringent

# Try MNN -----------------------------------------------------------------

jmark1 <- "H3K4me3"
# jmark2 <- "H3K9me3"
jmark2 <- "H3K27me3"
# jmark2 <- "H3K27me3"
# jmark2 <- "H3K4me3"

keeptop <- 100
# filter variable genes
top.genes.X <- apply(tm.result.lst[[jmark1]]$terms, 1, function(jrow){
  indx <- sort(jrow, decreasing = TRUE, index.return=TRUE)
  return(indx$ix[1:keeptop])
}) %>%
  as.data.frame() %>%
  unlist() %>%
  unique()

top.genes.Y <- apply(tm.result.lst[[jmark2]]$terms, 1, function(jrow){
  indx <- sort(jrow, decreasing = TRUE, index.return=TRUE)
  return(indx$ix[1:keeptop])
}) %>%
  as.data.frame() %>%
  unlist() %>%
  unique()
top.genes.common <- intersect(top.genes.X, top.genes.Y)

X <- t(tm.result.lst[[jmark1]]$topics %*% tm.result.lst[[jmark1]]$terms[, top.genes.common])
Y <- t(tm.result.lst[[jmark2]]$topics %*% tm.result.lst[[jmark2]]$terms[, top.genes.common])

# find neighbors that are farthest in umap1 left

# plot granu
jcell1 <- (subset(dat.umap.long.trajs[[jmark1]]) %>% dplyr::filter(umap2 == min(umap2)))$cell
jcell2<- (subset(dat.umap.long.trajs[[jmark2]]) %>% dplyr::filter(umap2 == max(umap2)))$cell

plot(X[, jcell1], Y[, jcell2])

panel <- list(X, Y)
names(panel) <- c(jmark1, jmark2)
panel.preprocessed <- lapply(panel, basicSeuratProc)
# panel.preprocessed <- lapply(panel, CreateSeuratObject)

conos.outf <- paste0("/Users/yeung/data/scchic/pdfs/integration/conos_integrate_", paste(jmark1, jmark2, sep = "-"), ".pdf")
if (!file.exists(conos.outf)){
  
  con <- Conos$new(panel.preprocessed, n.cores=4)
  str(con$samples,1)
  # con$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=6)
  # con$buildGraph(k=30, k.self=5, space='CCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=FALSE, verbose=TRUE)
  con$buildGraph(k=30, k.self=5, space='CCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=FALSE, verbose=TRUE)
  con$findCommunities(method=leiden.community, resolution=1)
  con$plotPanel(font.size=4)
  con$plotGraph(alpha=0.1)
  
  saveRDS(con, file = conos.outf)
} else {
  con <- readRDS(conos.outf)
}

# plot output
con.out <- data.frame(cell = rownames(con$embedding), embed1 = con$embedding[, 1], embed2 = con$embedding[, 2])
con.clusters <- data.frame(cell = names(con$clusters$leiden$groups), embed.groups = con$clusters$leiden$groups)
con.out <- left_join(con.out, con.clusters)
con.out <- left_join(con.out, bind_rows(dat.umap.long.trajs[[jmark1]], dat.umap.long.trajs[[jmark2]]))

ggplot(con.out, aes(x = embed1, y = embed2, color = eryth)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(con.out, aes(x = umap1, y = umap2, color = embed.groups)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~mark)

# 
# # Jointly  ----------------------------------------------------------------
# 

imput.lst <- lapply(jmarks, function(jmark){
  basicSeuratProc(t(tm.result.lst[[jmark]]$topics %*% tm.result.lst[[jmark]]$terms))
})

conos.all.outf <- paste0("/Users/yeung/data/scchic/pdfs/integration/conos_integrate_all.rds")
# conos.all.rds.outf <- paste0("/Users/yeung/data/scchic/pdfs/integration/conos_integrate_all.rds")
con.all <- Conos$new(imput.lst, n.cores=4)
str(con.all$samples, 1)
# con.all$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=6)
con.all$buildGraph(k=30, k.self=5, space='CCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=FALSE, verbose=TRUE)
con.all$findCommunities(method=leiden.community, resolution=1)
con.all$plotPanel(font.size=4)
con.all$plotGraph(alpha=0.1)
saveRDS(con.all, file = conos.all.outf)

con.out.all <- data.frame(cell = rownames(con.all$embedding), embed1 = con.all$embedding[, 1], embed2 = con.all$embedding[, 2])
con.out.all.clusters <- data.frame(cell = names(con.all$clusters$leiden$groups), embed.groups = con.all$clusters$leiden$groups)
con.out.all <- left_join(con.out.all, con.out.all.clusters)
con.out.all <- left_join(con.out.all, bind_rows(dat.umap.long.trajs))
# con.all.merge <- left_join(con.all, bind_rows(dat.umap.long.trajs))

ggplot(con.out.all, aes(x = umap1, y = umap2, color = embed.groups)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette) + 
  facet_wrap(~mark)



# MNN ---------------------------------------------------------------------

# 
# rds.path <- paste0("/Users/yeung/data/scchic/pdfs/integration/mnn_out.", paste(jmark1, jmark2, sep = "-"), ".rds")
# if (!file.exists(rds.path)){
#   system.time(
#     out <- scran::fastMNN(X, Y)
#   )
#   saveRDS(out, file=rds.path)
# } else {
#   print("file alrdy exists, loading it")
#   out <- readRDS(rds.path)
# }
# 
# plot(out$corrected[, 1], out$corrected[, 2])
# 
# # explore nearest neighbors
# mnn.cor <- data.frame(cell = c(colnames(X), colnames(Y)), PC1 = out$corrected[, 1], PC2 = out$corrected[, 2], indx = seq(nrow(out$corrected)))
# 
# # check if pairs are matching?
# 
# pair.cells.indx <- unlist(as.data.frame(out$pairs[[1]]))
# 
# mnn.corr.merged <- left_join(mnn.cor, bind_rows(dat.umap.long.trajs[[jmark1]], dat.umap.long.trajs[[jmark2]]))
# 
# mnn.corr.merged$is.mnn <- mnn.corr.merged$indx %in% pair.cells.indx
# 
# ggplot(mnn.corr.merged, aes(x = PC1, y = PC2, color = eryth)) + geom_point()  + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~mark)
# 
# ggplot(mnn.corr.merged, aes(x = PC1, y = PC2, color = is.mnn)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~mark)

# get neighbors
# 
# pdf(file = paste0("/Users/yeung/data/scchic/pdfs/integration/integrate_", jmark1, "-", jmark2, ".pdf"), useDingbats = FALSE)
# for (i in seq(nrow(out$pairs[[1]]))[1:200]){
#   print(i)
#   cells.indx.filt <- unlist(out$pairs[[1]][i, ])
#   jcells <- mnn.corr.merged[cells.indx.filt, ]$cell
#   m <- ggplot(mnn.corr.merged %>% mutate(is.pair = cell %in% jcells) %>% arrange(is.pair),
#               aes(x = umap1, y = umap2, color = is.pair, size = is.pair)) + geom_point()  + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     facet_wrap(~mark)
#   print(m)
# }
# dev.off()

# # what is relationship between two granu cells?
# jcell1 <- (subset(mnn.corr.merged, mark == jmark1 & is.mnn) %>% dplyr::filter(umap1 == max(umap1)))$cell
# jcell2 <- (subset(mnn.corr.merged, mark == jmark2 & is.mnn) %>% dplyr::filter(umap2 == max(umap2)))$cell
# 
# plot(unlist(X[, jcell1]), unlist(Y[, jcell2]))



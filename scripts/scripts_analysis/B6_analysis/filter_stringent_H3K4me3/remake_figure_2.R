# Jake Yeung
# Date of Creation: 2019-06-05
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/remake_figure_2.R
# Remake figures for main Figure 2.
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)
source("scripts/Rfunctions/PlotFunctions.R")



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- "H3K4me3"
# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.RData"
load(inf, v=T)

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)
out.objs.stringent <- out.objs
tm.result.stringent <- posterior(out.objs$out.lda)

inf.lda.all <- "/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
load(inf.lda.all, v=T)


out.objs$H3K4me3 <- out.objs.stringent

inf.objs <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.objs, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent$H3K4me3

# Plot the 4 umaps  -------------------------------------------------------

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K27me3 = "darkorange1", H3K9me3 = "red1")

ncells <- lapply(jmarks, function(jmark) return(length(unique(dat.umap.long.trajs[[jmark]]$cell))))

mlst <- lapply(jmarks, function(jmark) PlotXYNoColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", jcol = colsvec[[jmark]], jsize = 1) + ggtitle(ncells[[jmark]]))





# Plot hits ---------------------------------------------------------------

# topics.vec <- c(7, 14, 40, 47, 48)
topics.vec <- c(34, 42, 39, 14, 40, 27)
hits <- c("Hbb-bs", "Il2ra", "Prf1", "S100a8", "Kit", "Gzmb", "Inpp4b", "Irf8", "S100a7a", "S100a6", "Tal1", "Sox6", "Hbb-y", "Gypa", "Bach2", "Sirpa", "Atf3", "Ccl2")
topname <- c("Erythroblasts", "B cells", "NK cells", "Granulocytes", "Progenitors")

topics.mat <- as.data.frame(posterior(out.objs[[jmark]]$out.lda)$topics)
colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))

# annotate distance to terms.filt
terms.filt.dist <- left_join(terms.filt, out2.df.closest %>% dplyr::select(region_coord, dist.to.tss, gene) %>% dplyr::rename(term = region_coord)) %>%
  group_by(gene) %>%
  filter(abs(dist.to.tss) == min(abs(dist.to.tss))) %>%
  # filter(weight == max(weight)) %>%
  filter(rnk == min(rnk))

terms.keep <- (terms.filt.dist %>%
  filter(gene %in% hits))$term

imput.mat <- t(posterior(out.objs[[jmark]]$out.lda)$topics %*% posterior(out.objs[[jmark]]$out.lda)$terms[, terms.keep])

imput.sub <- tidyr::gather(data.frame(term = terms.keep, imput.mat[terms.keep, ]), key = "cell", value = "exprs", -term)
imput.sub$cell <- gsub("\\.", "-", imput.sub$cell)
imput.sub <- left_join(imput.sub, subset(terms.filt.dist, term %in% terms.keep))

jscale <- 10^6
jpseudo <- 0
do.log <- FALSE
if (do.log){
  imput.sub <- imput.sub %>%
    mutate(exprs = log2(exprs * jscale + jpseudo)) 
}

imput.wide <- tidyr::spread(imput.sub %>% dplyr::select(gene, exprs, cell), key = "gene", value = "exprs")

topics.mat$cell <- rownames(topics.mat)

# for topics and gene
dat.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat %>% dplyr::select(c("cell", paste0("topic_", topics.vec))))
dat.merge <- left_join(dat.merge, imput.wide)

jsub <- dat.merge

steps <- c("lightblue", "gray95", scales::muted("red"))
steps2 <- c("blue", "white", "red")
steps.top <- c("gray95", "gray50", scales::muted("blue"))
pal <- color.palette(steps, c(10, 10), space="rgb")
pal2 <- color.palette(steps2, c(10, 10), space="rgb")
pal.top <- color.palette(steps.top, c(10, 10), space="rgb")

jsize <- 2

pdf("/Users/yeung/data/scchic/pdfs/B6_figures/stringent/redo_figure2.pdf", useDingbats = FALSE)
multiplot(mlst[[1]], mlst[[3]], mlst[[2]], mlst[[4]], cols = 2)
for (hit in hits){
  print(hit)
  jtitle <- hit
  jcname <- hit
  jsub <- RankOrder(jsub, cname = jcname, out.cname = "orderrank")
  jsub.log <- jsub
  jsub.log[[hit]] <- log2(jsub[[hit]] * jscale + jpseudo)
  xvar <- "umap1"; yvar <- "umap2"; cname.str <- jcname
  m1 <- ggplot(jsub, aes_string(x = xvar, y = yvar, col = paste0("`", jcname, "`"), order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = pal(50))
  print(m1)
  m1.log <- ggplot(jsub.log, aes_string(x = xvar, y = yvar, col = paste0("`", jcname, "`"), order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = pal(50))
  print(m1.log)
  m1.br <- ggplot(jsub, aes_string(x = xvar, y = yvar, col = paste0("`", jcname, "`"), order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = pal2(50))
  print(m1.br)
  m1.br.log <- ggplot(jsub.log, aes_string(x = xvar, y = yvar, col = paste0("`", jcname, "`"), order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = pal2(50))
  print(m1.br.log)
}

# plot topics
for (jtop in topics.vec){
  topicname <- paste0("topic_", jtop)
  jtitle <- topicname
  jcname <- topicname
  jsub <- RankOrder(jsub, cname = jcname, out.cname = "orderrank")
  xvar <- "umap1"; yvar <- "umap2"; cname.str <- jcname
  m1.top <- ggplot(jsub, aes_string(x = xvar, y = yvar, col = paste0("`", jcname, "`"), order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = pal.top(50))
  print(m1.top)
}

# print umap with clusters
jsize.umap <- 2.5
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
m1.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point(size = jsize.umap) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
print(m1.louvain)
m1.louvain2 <- PlotXYWithColor(dat.umap.long, xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jsize = jsize.umap)
print(m1.louvain2)

dev.off()



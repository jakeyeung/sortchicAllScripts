# Jake Yeung
# Date of Creation: 2019-05-21
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig2.R
# Load objects and make fig 2


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)

source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")

# Constants ---------------------------------------------------------------

datmain <- "/Users/yeung/data/scchic/robjs/B6_objs"
outdirplot <- "/Users/yeung/data/scchic/pdfs/B6_figures/trajectories"
inmain <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"

# add trajectory
inf.traj <- file.path(datmain, "traj_objs_all_marks.Rdata")
assertthat::assert_that(file.exists(inf.traj))

inf.objs <- file.path(datmain, "LDA_objects_all_marks.Rdata")
assertthat::assert_that(file.exists(inf.objs))

inf.termsfilt <- file.path(datmain, "terms_filt_H3K4me1_bin_TRUE_k_50.genomewide.RData")
assertthat::assert_that(file.exists(inf.termsfilt))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- jmarks[["H3K4me1"]]
inf.dats <- lapply(jmarks, function(jmark) file.path(inmain, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))



lapply(inf.dats, function(x) assertthat::assert_that(file.exists(x)))

topics.vec <- c(7, 14, 40, 47, 48)
hits <- c("Hbb-bs", "Il2ra", "Prf1", "S100a8", "Kit", "Gzmb", "Inpp4b", "Irf8", "S100a7a", "S100a6", "Tal1", "Sox6")
topname <- c("Erythroblasts", "B cells", "NK cells", "Granulocytes", "Progenitors")




# Load umaps  -------------------------------------------------------------

load(inf.traj, v=T)

load(inf.termsfilt, v=T)  # terms.filt

load(inf.objs, v=T)

dat.umap.long <- lapply(inf.dats, LoadUmap) %>%
  bind_rows() %>%
  dplyr::select(-repl, -techname)


# Make 4 heatmaps different colors each  ----------------------------------

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K27me3 = "darkorange1", H3K9me3 = "red1")

mlst <- lapply(jmarks, function(jmark) PlotXYNoColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", jcol = colsvec[[jmark]], jsize = 1))

pdf(file = "/tmp/fig_2_output.pdf", useDingbats = FALSE)

multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)


# Make topic loadings -----------------------------------------------------

# H3K4me1 only
# show topic loadings
jmark <- "H3K4me1"
topics.mat <- data.frame(cell = rownames(out.objs[[jmark]]$tm.result$topics), out.objs[[jmark]]$tm.result$topics, stringsAsFactors = FALSE)

dat.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat)

# replace X with Topic_ in colnames
colnames(dat.merge) <- gsub("^X", "Topic_", colnames(dat.merge))

# select from https://stats.biopapyrus.jp/r/graph/rcolorbrewer.html
# or https://stats.biopapyrus.jp/r/graph/rcolorbrewer.html
xvar <- "umap1"; yvar <- "umap2"
jsize <- 3
for (jtopic in topics.vec){
  jsub <- dat.merge
  jtitle <- jtopic
  jcname <- paste0("Topic_", jtopic)
  cname.str <- paste0("Topic_", jtopic)
  jsub <- RankOrder(jsub, cname = jcname, out.cname = "orderrank")
  m1 <- ggplot(jsub, aes_string(x = xvar, y = yvar, col = jcname, order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(7, "Reds"))
    # scale_color_gradientn(colours = RColorBrewer::brewer.pal(7, "YlOrRd"))
    # scale_color_gradient(low = "gray85", high = scales::muted("darkred"))
  print(m1)
}


# eryth, bcells, nk cells, granu, progenitors
topics.vec <- c(7, 14, 40, 47, 48)
hits <- c("Hbb-bs", "Il2ra", "Prf1", "S100a8", "Kit", "Gzmb", "Inpp4b", "Irf8", "S100a7a", "S100a6", "Tal1", "Sox6", "Hbb-y", "Gypa")
topname <- c("Erythroblasts", "B cells", "NK cells", "Granulocytes", "Progenitors")


# Make gene examples ------------------------------------------------------

# H3K4me1 only 
topics.mat <- as.data.frame(posterior(out.objs[[jmark]]$out.lda)$topics)
colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))

imput.mat <- t(posterior(out.objs[[jmark]]$out.lda)$topics %*% posterior(out.objs[[jmark]]$out.lda)$terms)

terms.keep <- subset(terms.filt, gene %in% hits)$term


imput.sub <- tidyr::gather(data.frame(term = terms.keep, imput.mat[terms.keep, ]), key = "cell", value = "exprs", -term)
imput.sub$cell <- gsub("\\.", "-", imput.sub$cell)
imput.sub <- left_join(imput.sub, terms.filt)
imput.wide <- tidyr::spread(imput.sub %>% dplyr::select(gene, exprs, cell), key = "gene", value = "exprs")

topics.mat$cell <- rownames(topics.mat)

# for topics and gene
dat.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat %>% dplyr::select(c("cell", paste0("topic_", topics.vec))))
dat.merge <- left_join(dat.merge, imput.wide)

jsub <- dat.merge
for (hit in hits){
  jtitle <- hit
  jcname <- hit
  jsub <- RankOrder(jsub, cname = jcname, out.cname = "orderrank")
  xvar <- "umap1"; yvar <- "umap2"; cname.str <- jcname
  m1 <- ggplot(jsub, aes_string(x = xvar, y = yvar, col = paste0("`", jcname, "`"), order = "orderrank")) + 
    ggrastr::geom_point_rast(size = jsize) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  + 
    xlab("") + ylab("") + ggtitle(jtitle) + 
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(7, "Blues"))
  # scale_color_gradientn(colours = RColorBrewer::brewer.pal(7, "YlOrRd"))
  # scale_color_gradient(low = "gray85", high = scales::muted("darkred"))
  print(m1)
}


dev.off()




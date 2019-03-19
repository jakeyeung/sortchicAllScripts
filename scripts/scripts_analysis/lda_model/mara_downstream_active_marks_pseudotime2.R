# Jake Yeung
# Date of Creation: 2019-03-07
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model/mara_downstream_active_marks_pseudotime2.R
# Do pseudotime on H3K4me1

rm(list=ls())

setwd("~/projects/scchic")
jtime <- Sys.time() 

library(data.table)
library(tidyr)
library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(forcats)
library(ggrepel)
library(biomaRt)

library(igraph)  # louvain

library(Gviz)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

source("scripts/Rfunctions/MaraDownstream.R")

source("scripts/Rfunctions/PlotFunctions.R")

GetTopGenes <- function(dat.long, jcelltype, jtopic, topn = 20){
  jsub <- dat.long %>% 
    rowwise() %>%
    mutate(CellType2 = ifelse(CellType == "nucleate_erythrocyte", TRUE, FALSE)) %>% 
    group_by(CellType2, Gene_Name) %>%
    summarise(zscore = mean(logFPKM)) %>%
    group_by(Gene_Name) %>% 
    summarise(zscore.diff = diff(range(zscore)))
  # jsub.hash <- hash(jsub$Gene_Name, jsub$zscore.diff)
  # get top genes
  jsub.top <- top.peaks.annotated %>% 
    filter(topic == jtopic & rnk <= 100) %>%
    group_by(topic, SYMBOL) %>%
    filter(beta == max(beta)) %>%
    arrange(rnk)
  jsub.top <- left_join(jsub.top, jsub, by = c("SYMBOL" = "Gene_Name")) %>% 
    arrange(desc(zscore.diff))
  jsub.sub <- jsub.top[1:topn, ]
  # genes <- jsub.top$SYMBOL[1:topn]
  return(jsub.sub)
}

# do Principal Curves on filtered genes
InferTrajOnUmap <- function(dat, cname = "granu", init.on = "umap1"){
  inmat <- subset(dat[which(dat[[cname]] == TRUE), ], select = c(umap1, umap2, cell))
  rownames(inmat) <- inmat$cell
  inmat$cell <- NULL
  # order from origin outwards
  inmat <- inmat[order(inmat$umap1 ^ 2 + inmat$umap2 ^ 2), ]
  # create initial ordering
  if (init.on == "umap1"){
    start.x <- inmat[["umap1"]]
    start.y <- rep(mean(inmat[["umap2"]]), nrow(inmat))
  } else if (init.on == "umap2"){
    start.x <- rep(mean(inmat[["umap1"]], nrow(inmat)))
    start.y <- inmat[["umap2"]]
  } else {
    warning(paste("init on must be umap1 or umap2", init.on))
  }
  startvals <- as.matrix(data.frame(umap1 = start.x, umap2 = start.y))
  proj <- principal_curve(as.matrix(inmat))
  dat.pca.proj <- as.data.frame(proj$s)
  dat.pca.proj$cell <- rownames(proj$s)
  # scale lambda from 0 to 1
  dat.pca.proj$lambda <- proj$lambda / max(proj$lambda)
  dat.pca.proj <- dat.pca.proj %>% arrange(lambda)
  # if end of lambda is near origin, then switch
  # endpoint <- dat.pca.proj[, 1]$umap1 *
  return(dat.pca.proj)
}

# Load bulk RNAseq data ---------------------------------------------------

dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))


# Load MARA activities  ---------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3")

mdirs <- lapply(jmarks, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.repress <- c("H3K27me3", "H3K9me3")
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  print(mdir)
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all

mdirs <- c(mdirs, mdirs.repress)
mara.outs <- lapply(mdirs, LoadMARA)

head(mara.outs[[1]]$act.long)


act.long.merged <- rbind(mara.outs[[1]]$act.long, mara.outs[[2]]$act.long, mara.outs[[3]]$act.long, mara.outs[[4]]$act.long)
# zscores.merged <- left_join(mara.outs[[3]]$zscores, mara.outs[[4]]$zscores, by = "motif")
# zscores.merged <- left_join(mara.outs[[3]]$zscores, mara.outs[[4]]$zscores, by = "motif")
zscores.merged <- purrr::reduce(list(mara.outs[[1]]$zscores, mara.outs[[2]]$zscores, mara.outs[[3]]$zscores, mara.outs[[4]]$zscores), left_join, by = "motif")
cnames <- c("motif", paste("zscore", jmarks.all, sep = "."))

colnames(zscores.merged) <- cnames


zscore.thres <- 0.75
zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
  ifelse(max(jrow[[2]], jrow[[3]]) > zscore.thres, jrow[[1]], NA)
})

zscores.merged.mean <- zscores.merged
zscores.merged.mean$zscore.mean <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, mean)
zscores.merged.mean$zscore.max <- apply(zscores.merged[, paste("zscore", jmarks.all, sep = ".")], 1, max)

zscores.merged.mean <- zscores.merged.mean %>%
  arrange(desc(zscore.mean))


# Constants you can tweak -------------------------------------------------

top.thres <- 0.995
inmain <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"
meanfilt <- 10
jbin <- "TRUE"; kstr <- "15_20_25_30_35"
# jbin <- "FALSE"; kstr <- "15_20_25_30"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
jmark <- "H3K4me1"
# jmark <- "H3K9me3"
indir <- paste0("lda_outputs.meanfilt_", 
                meanfilt,
                ".cellmin_100.cellmax_500000.binarize.", jbin, 
                ".no_filt")
fname <- paste0("lda_out_meanfilt.BM-", jmark, ".CountThres0.K-", kstr, ".Robj")
inf <- file.path(inmain, indir, fname)
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

out.lda <- ChooseBestLDA(out.lda)
(kchoose <- out.lda@k)
tm.result <- posterior(out.lda)

outmain <- "~/Dropbox/scCHiC_figs/FIG4_BM/primetime_plots"
plotout <- file.path(outmain, paste0(jmark, "_LDA_summarized_bins.pdf"))

jtopic <- 1


topics.mat <- tm.result$topics
terms.mat <- tm.result$terms

nn=40
nnterms <- 15
jmetric='euclidean' 
jmindist=0.4
jseed=123

custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed)
custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# check your umap settings
jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)

print(sessionInfo())

# Plot dat umap -----------------------------------------------------------
# jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("lightblue", "darkblue"))
jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, topics.mat, c("pink", "darkred"))
nb.col <- 5
nb.row <- ceiling(kchoose / nb.col)

plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, asp = 0.75, cex = 0.2)
par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1), pty = "s")
mapply(function(jcol.rgb, jtopic){
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.75, cex = 0.2)
}, jcol.rgbs, seq(kchoose))

# Plot terms umap ---------------------------------------------------------
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
top.regions <- unique(unlist(topic.regions))
terms.mat <- t(tm.result$terms)[top.regions, ]

# Uncomment below to plot the UMAP on the terms, can take a few more minutes

# dat.umap.terms <- umap(terms.mat, config = custom.settings.terms)
# # downsample rows for plotting purposes
# downsamp.i <- sample(seq(nrow(dat.umap.terms$layout)), size = round(0.1 * nrow(dat.umap.terms$layout)), replace = FALSE)
# jcol.rgb.terms <- lapply(seq(kchoose), ColorsByGamma, terms.mat[downsamp.i, ], c("pink", "red", "darkred"))
# par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
# mapply(function(jcol.rgb.term, jtopic){
#   plot(dat.umap.terms$layout[downsamp.i, 1], dat.umap.terms$layout[downsamp.i, 2],
#        col = jcol.rgb.term, pch = 20, asp = 0.75,
#        main = paste("Peak Weights, T", jtopic))
# }, jcol.rgb.terms, seq(kchoose))

# how many cells in the island?
cell.assign <- apply(topics.mat, 1, which.max)

top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>% 
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta)))

# annotate regions?
regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd), 
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)

top.peaks.annotated <- dplyr::left_join(top.peaks, subset(regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))

hit.peaks <- subset(regions.annotated, abs(distanceToTSS) < 50000 & grepl("Hbb", SYMBOL))$region_coord

jpeak <- "chr7:103800000-103900000"
PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = TRUE, return.plot.only = TRUE, usettings=custom.settings)

# how top hits for specific topics
# Progenitor cellsare in Topic 12?
topics.mat.named <- as.data.frame(topics.mat)
topics.mat.named$cell <- rownames(topics.mat.named)

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout))
dat.umap.long <- left_join(dat.umap.long, topics.mat.named)

# plot northern island
jcol.rgb <- jcol.rgbs[[jtopic]]
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# make it pretty

m.all <- ggplot(dat.umap.long, aes_string(x = "umap1", y = "umap2")) + geom_point(size = 0.8, alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_continuous(paste0("Topic ", jtopic, "\nWeight"))
print(m.all)

# plot top hits for topic 12
print(subset(top.peaks.annotated, topic == 12), n = 50)
print(subset(top.peaks.annotated, topic == 10), n = 50)
print(subset(top.peaks.annotated, topic == 7), n = 50)

jgene <- c("Tgm2", "Mmrn1", "Trim48", "Pdzk1ip1", "Mllt3", "Mecom")
jgene <- c("Vldlr", "Uhrf1", "Rrm2", "Lig1", "Tipin")
jgene <- c("F2r", "Itga2b", "Zfp385a", "Zfpm1", "Plek", "Cd9", "Zeb2")

# translate beta to log fold change?
mat.norm <- t(tm.result$topics %*% tm.result$terms)  # this should give normalized signal, without the poisson noise?
# mat.norm <- mat.norm[top.regions, ]

# label using louvain clustering?

# do new knn

# dat.umap.fo

nn.louv=27
# nn.louv=15
jmetric.louv='euclidean' 
jmindist.louv=0.4
jseed.louv=123

custom.settings.louv <- GetUmapSettings(nn=nn.louv, 
                                        jmetric=jmetric.louv, 
                                        jmindist=jmindist.louv, 
                                        seed = jseed.louv)

dat.umap.louv <- umap(topics.mat, config = custom.settings.louv)

dat.umap.louv.long <- data.frame(umap1 = dat.umap.louv$layout[, 1], umap2 = dat.umap.louv$layout[, 2], cell = rownames(dat.umap.louv$layout), 
                                 stringsAsFactors = FALSE)

cell.indx <- hash(rownames(dat.umap.louv$knn$indexes), dat.umap.louv$knn$indexes[, 1])
cell.indx.rev <- hash(dat.umap.louv$knn$indexes[, 1], rownames(dat.umap.louv$knn$indexes))
nr <- nrow(dat.umap.louv$knn$indexes)
nc <- ncol(dat.umap.louv$knn$indexes)
edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
colnames(edgelist) <- c("from", "to")
for (vertex.i in seq(nr)){
  istart <- nc*(vertex.i - 1)+1
  iend <- nc*vertex.i
  edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
  edgelist[istart : iend, 2] <- sapply(dat.umap.louv$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
  # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
}
g <- graph_from_data_frame(edgelist, directed=FALSE)
g.out <- cluster_louvain(g, weights = NULL)
V(g)$color <- g.out$membership
clstr <- hash(g.out$names, g.out$membership)

dat.umap.long$louvain <- sapply(dat.umap.long$cell, function(x) clstr[[x]])

jclst <- 6
clstrs.orig <- as.character(sort(unique(as.numeric(dat.umap.long$louvain))))
# swap jclst with first element
clstrs.new <- clstrs.orig
clstrs.new[c(1, which(clstrs.new == jclst))] <- clstrs.new[c(which(clstrs.new == jclst), 1)]
remap.clstr <- hash(clstrs.orig, clstrs.new)


dat.umap.long$louvain <- sapply(as.character(dat.umap.long$louvain), function(x) remap.clstr[[x]])
dat.umap.long$louvain <- factor(as.character(dat.umap.long$louvain), levels = clstrs.orig)  # 1 to N
m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Spectral")
print(m.louvain)

# plot graph with edges?
# https://stackoverflow.com/questions/5364264/how-to-control-the-igraph-plot-layout-with-fixed-positions

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
coords <- layout.auto(g)
plot.igraph(igraph::simplify(g),
            layout = dat.umap$layout[V(g)$name, ],
            vertex.label = NA,
            edge.curved=FALSE,
            label = NA,
            edge.width = 0.5,
            vertex.size = 1)


# Plot example genes ------------------------------------------------------

jsize <- 1
jcolvec <- c("blue", "gray80", "red")

jpeak <- c("chr7:103820000-103920000")
jgene <- "Hbb-bs"

jgene <- "Cebpb"
(jpeak <- subset(top.peaks.annotated, grepl(jgene, SYMBOL))$term[[1]])
# jpeak <- "chr7:35080000-35180000"

m1 <- PlotImputedPeaks2(tm.result, jpeak, jmark,
                        use.count.mat = NULL,
                        usettings=dat.umap, 
                        gname = jgene,
                        jsize = jsize, jcolvec = jcolvec, .log = TRUE)
print(m1)

# project granulocyte cells onto pseudotime

print(m.louvain)

# show granulocyte trajectory
print(m.louvain)
# clsts <- c("9", "3", "7")
clsts <- c("7")
# clsts <- c("8", "2", "6")

m.louvain.filt <- ggplot(dat.umap.long %>% mutate(louvain = ifelse(louvain %in% clsts, TRUE, FALSE)), aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.louvain.filt)

# keep cells in granulocyte trajectory. Split louvain clstr 7 into two parts
dat.sub.granu <- dat.umap.long %>% filter((louvain == 7 & umap2 > 0) | (louvain %in% c("9", "3")))
cells.granu <- hash(dat.sub.granu$cell, TRUE)
dat.umap.long$granu <- sapply(as.character(dat.umap.long$cell), function(x) ifelse(!is.null(cells.granu[[x]]), cells.granu[[x]], FALSE))

# do erythryocytes easy
dat.sub.eryth <- dat.umap.long %>% filter(louvain == 1)
cells.eryth <- hash(dat.sub.eryth$cell, TRUE)
dat.umap.long$eryth <- sapply(as.character(dat.umap.long$cell), function(x) ifelse(!is.null(cells.eryth[[x]]), cells.eryth[[x]], FALSE))

# megakaryocytes
dat.sub.mega <- dat.umap.long %>% filter(louvain == 7)
cells.mega <- hash(dat.sub.mega$cell, TRUE)
dat.umap.long$mega <- sapply(as.character(dat.umap.long$cell), function(x) ifelse(!is.null(cells.mega[[x]]), cells.mega[[x]], FALSE))

# b lymphos
dat.sub.bcell <- dat.umap.long %>% filter(louvain %in% c("8", "2", "6"))
cells.bcell <- hash(dat.sub.bcell$cell, TRUE)
dat.umap.long$bcell <- sapply(as.character(dat.umap.long$cell), function(x) ifelse(!is.null(cells.bcell[[x]]), cells.bcell[[x]], FALSE))

# tcells split cluster 8 into two parts
dat.sub.tcell <- dat.umap.long %>% filter((louvain == 8 & umap1 > -3) | (louvain == 5))
cells.tcell <- hash(dat.sub.tcell$cell, TRUE)
dat.umap.long$tcell <- sapply(as.character(dat.umap.long$cell), function(x) ifelse(!is.null(cells.tcell[[x]]), cells.tcell[[x]], FALSE))


m.louvain.filt <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = granu)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.louvain.filt)


traj.granu <- InferTrajOnUmap(dat.umap.long, cname = "granu", init.on = "umap1")
traj.eryth <- InferTrajOnUmap(dat.umap.long, cname = "eryth", init.on = "umap2")
traj.bcell <- InferTrajOnUmap(dat.umap.long, cname = "bcell", init.on = "umap1")
traj.tcell <- InferTrajOnUmap(dat.umap.long, cname = "tcell", init.on = "umap1")
traj.mega <- InferTrajOnUmap(dat.umap.long, cname = "mega", init.on = "umap2")

# plot original
jsize <- 2
ggplot(dat.umap.long %>% mutate(louvain = ifelse(louvain %in% clsts, TRUE, FALSE)), aes(x = umap1, y = umap2)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = traj.granu, aes(color = lambda), size = jsize) + 
  geom_path(data = traj.eryth, aes(color = lambda), size = jsize) + 
  geom_path(data = traj.bcell, aes(color = lambda), size = jsize) + 
  geom_path(data = traj.tcell, aes(color = lambda), size = jsize) + 
  geom_path(data = traj.mega, aes(color = lambda), size = jsize) 


# Analyze on cell activities on pseudotime  -------------------------------


jmotif <- "Cebpb"
jgene <- "Cebpb"
(jpeak <- subset(top.peaks.annotated, grepl(jgene, SYMBOL))$term[[1]])

mat.imput <- t(tm.result$topics %*% tm.result$terms)
peak.counts <- data.frame(peak = jpeak, 
                          cell = colnames(mat.imput),
                          counts = log10(mat.imput[jpeak, ]))
act.counts <- subset(act.long.merged, motif == jmotif)

dat.counts <- left_join(peak.counts, act.counts)

act.counts$counts <- act.counts$activity

act.counts <- left_join(act.counts, traj.granu)
peak.counts <- left_join(peak.counts, traj.granu)

traj.granu2 <- traj.granu
traj.granu2 <- left_join(traj.granu2, dat.counts)

# show umap, ad cebpb activity 
dat.umap.long.cat <- left_join(dat.umap.long, act.counts %>% dplyr::select(cell, counts))
jsize <- 2
ggplot(dat.umap.long.cat, aes(x = umap1, y = umap2, color = counts)) + 
  geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_path(data = traj.granu, aes(x = umap1, y = umap2), size = jsize, alpha = 0.5, inherit.aes = FALSE)

# plot pseudotime
ggplot(traj.granu2, aes(x = lambda, y = counts, color = counts)) + geom_smooth(se = FALSE) + geom_point()
ggplot(traj.granu2, aes(x = lambda, y = activity, color = activity)) + geom_smooth(se = FALSE) + geom_point()

ggplot(act.counts %>% filter(!is.na(lambda)), aes(x = lambda, y = activity)) + geom_point() + geom_smooth(se = FALSE)
ggplot(peak.counts %>% filter(!is.na(lambda)), aes(x = lambda, y = counts)) + geom_point() + geom_smooth(se = FALSE)

# merge the two data and replot
dat.merge <- bind_rows(act.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "activity"), peak.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "counts"))
dat.merge <- dat.merge %>%
  group_by(assay) %>%
  mutate(counts = scale(counts))

ggplot(dat.merge, aes(x = lambda, y = counts, color = assay, group = assay)) + geom_point() + geom_smooth(se = FALSE)

# Analyze on cell activities on pseudotime: eryth  -------------------------------

jmotif <- "Tal1"
jgene <- "Tal1"

# jmotif <- "Pml"
# jgene <- "Pml"

(jpeak <- subset(top.peaks.annotated, grepl(jgene, SYMBOL))$term[[1]])

# mat.imput <- t(tm.result$topics %*% tm.result$terms)
peak.counts <- data.frame(peak = jpeak, 
                          cell = colnames(mat.imput),
                          counts = log10(mat.imput[jpeak, ]))
act.counts <- subset(act.long.merged, motif == jmotif)

dat.counts <- left_join(peak.counts, act.counts)

act.counts$counts <- act.counts$activity

act.counts <- left_join(act.counts, traj.eryth)
peak.counts <- left_join(peak.counts, traj.eryth)

# merge the two data and replot
dat.merge <- bind_rows(act.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "activity"), peak.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "counts"))
dat.merge <- dat.merge %>%
  group_by(assay) %>%
  mutate(counts = scale(counts))

ggplot(dat.merge, aes(x = lambda, y = counts, color = assay, group = assay)) + geom_point() + geom_smooth(se = FALSE)


# Cell activities B cells  ------------------------------------------------

jmotif <- "Nfatc1"
jgene <- "Nfatc3"

(jpeak <- subset(top.peaks.annotated, grepl(jgene, SYMBOL))$term[[1]])

# mat.imput <- t(tm.result$topics %*% tm.result$terms)
peak.counts <- data.frame(peak = jpeak, 
                          cell = colnames(mat.imput),
                          counts = log10(mat.imput[jpeak, ]))
act.counts <- subset(act.long.merged, motif == jmotif)

dat.counts <- left_join(peak.counts, act.counts)

act.counts$counts <- act.counts$activity

act.counts <- left_join(act.counts, traj.bcell)
peak.counts <- left_join(peak.counts, traj.bcell)

# merge the two data and replot
dat.merge <- bind_rows(act.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "activity"), peak.counts %>% filter(!is.na(lambda)) %>% mutate(assay = "counts"))
dat.merge <- dat.merge %>%
  group_by(assay) %>%
  mutate(counts = scale(counts))

ggplot(dat.merge, aes(x = lambda, y = counts, color = assay, group = assay)) + geom_point() + geom_smooth(se = FALSE)


# Make plots --------------------------------------------------------------

# make through louvain 

nn.louv=65
jmetric.louv='euclidean' 
jmindist.louv=0.2
jseed.louv=123

custom.settings.louv <- GetUmapSettings(nn=nn.louv, 
                                        jmetric=jmetric.louv, 
                                        jmindist=jmindist.louv, 
                                        seed = jseed.louv)

dat.umap.louv2 <- umap(topics.mat, config = custom.settings.louv)

dat.umap.louv.long2 <- data.frame(umap1 = dat.umap.louv2$layout[, 1], umap2 = dat.umap.louv2$layout[, 2], cell = rownames(dat.umap.louv2$layout), 
                                 stringsAsFactors = FALSE)

cell.indx <- hash(rownames(dat.umap.louv2$knn$indexes), dat.umap.louv2$knn$indexes[, 1])
cell.indx.rev <- hash(dat.umap.louv2$knn$indexes[, 1], rownames(dat.umap.louv2$knn$indexes))
nr <- nrow(dat.umap.louv2$knn$indexes)
nc <- ncol(dat.umap.louv2$knn$indexes)
edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
colnames(edgelist) <- c("from", "to")
for (vertex.i in seq(nr)){
  istart <- nc*(vertex.i - 1)+1
  iend <- nc*vertex.i
  edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
  edgelist[istart : iend, 2] <- sapply(dat.umap.louv$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
  # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
}
g <- graph_from_data_frame(edgelist, directed=FALSE)
g.out <- cluster_louvain(g, weights = NULL)
V(g)$color <- g.out$membership
clstr <- hash(g.out$names, g.out$membership)

dat.umap.long2 <- dat.umap.long
dat.umap.long2$louvain <- sapply(dat.umap.long2$cell, function(x) clstr[[x]])

jclst <- 6
clstrs.orig <- as.character(sort(unique(as.numeric(dat.umap.long2$louvain))))
# swap jclst with first element
clstrs.new <- clstrs.orig
# clstrs.new <- clstrs.orig
# clstrs.new[c(1, which(clstrs.new == jclst))] <- clstrs.new[c(which(clstrs.new == jclst), 1)]
remap.clstr <- hash(clstrs.orig, clstrs.new)


dat.umap.long2$louvain <- sapply(as.character(dat.umap.long2$louvain), function(x) remap.clstr[[x]])
dat.umap.long2$louvain <- factor(as.character(dat.umap.long2$louvain), levels = clstrs.orig)  # 1 to N
m.louvain2 <- ggplot(dat.umap.long2, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_brewer(palette = "Spectral")
print(m.louvain2)



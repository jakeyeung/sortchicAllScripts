# Jake Yeung
# Date of Creation: 2019-05-29
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/10-make_primetime_figures_H3K4me3.R
# Do it for H3K4me3

rm(list=ls())


library(JFuncs)
library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)
library(ggrepel)

library(tidyr)

library(hash)
library(igraph)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Constants ---------------------------------------------------------------


jmark <- "H3K4me3"
keep.top.genes <- 150
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue",  H3K27me3 = "darkorange1", H3K9me3 = "red1")

# jbin <- TRUE; kstr <- "25_30_40_50"
jbin.vec <- c(TRUE, TRUE, FALSE, FALSE)
kstr.vec <- c("25_30_40_50", "25_30_40_50", "30_40_50", "30_40_50")
names(jbin.vec) <- jmarks
names(kstr.vec) <- jmarks

# jmark <- "H3K4me1"

# Load term sfilt ---------------------------------------------------------


load(paste0("/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_", jmark, "_bin_TRUE_k_50.genomewide_nofilt.RData"), v=T)  # terms.filt
terms.all <- terms.filt

load(paste0("/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_", jmark, "_bin_TRUE_k_50.genomewide.RData"), v=T)  # terms.filt
# terms.all <- terms.filt

# Load trajs --------------------------------------------------------------

load("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata", v=T)


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load files --------------------------------------------------------------




infs <- mapply(function(jbin, jmark, kstr) paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.", jbin, ".no_filt/lda_out_meanfilt.B6_", jmark, "_pcutoff_0.CountThres0.K-", kstr, ".Robj"), 
                  jbin.vec, jmarks, kstr.vec)
lapply(infs, function(inf) assertthat::assert_that(file.exists(inf)))


kchoose <- "50"
out.objs <- lapply(jmarks, function(jmark) LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = infs[[jmark]], convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose))
# out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose), jmarks, infs, SIMPLIFY = FALSE)
# names(out.objs) <- jmarks
# print(paste("K:", out.objs$out.lda@k))


# Add zscore info for cell typing -----------------------------------------
# 
# tm.result <- out.objs[[jmark]]$tm.result
# # get top topics 
# # analyze topic matrix across cells
# topics.mat <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics))
# topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>%
#   gather(key = "topic", value = "weight", -cell) %>%
#   rowwise() %>%
#   mutate(topic = as.numeric(substr(topic, 2, nchar(topic)))) %>%
#   group_by(topic) %>%
#   mutate(zscore = scale(weight, center = TRUE, scale = TRUE))
# 
# topics.sum <- topics.long %>%
#   group_by(topic) %>% # do entropy on 1 to 99% of cells
#   filter(zscore < quantile(zscore, 0.97)) %>%
#   mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
#   summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
#   arrange(entropy)
# print(topics.sum)
# 
# terms.long <- data.frame(term = colnames(tm.result$terms), as.data.frame(t(tm.result$terms)), stringsAsFactors = FALSE) %>%
#   gather(key = "topic", value = "weight", -term) %>%
#   mutate(topic = gsub("X", "", topic)) %>%
#   group_by(topic) %>%
#   arrange(desc(weight)) %>%
#   mutate(rnk = seq(length(weight))) %>%
#   rowwise()


tm.result <- out.objs[[jmark]]$tm.result
topics.long <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics)) %>%
  gather(key = "topic", value = "weight", -cell) %>%
  rowwise() %>%
  mutate(topic = as.numeric(substr(topic, 2, nchar(topic)))) %>%
  group_by(topic) %>%
  mutate(zscore = scale(weight, center = TRUE, scale = TRUE))

topics.sum <- topics.long %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(zscore < quantile(zscore, 0.97)) %>%
  mutate(zscore.prob = exp(zscore) / sum(exp(zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)
print(topics.sum)

topics.all <- topics.sum$topic
names(topics.all) <- topics.all

top.genes.all.lst <- lapply(topics.all, function(jtopic){
  jtmp <- subset(terms.filt, topic == jtopic)$gene
  jtmp <- jtmp[1:min(length(jtmp), keep.top.genes)]
  return(jtmp)
})

# get averag zscore across celltypes for each topic
zscores.avg.all <- lapply(topics.all, function(jtopic){
  jtopic <- as.character(jtopic)
  top.genes <- top.genes.all.lst[[jtopic]]
  jsub <- subset(dat.long, Gene_Name %in% top.genes) %>%
    group_by(CellType) %>%
    summarise(zscore = mean(zscore)) %>%
    ungroup() %>%
    mutate(weight = exp(zscore) / sum(exp(zscore)))
  jsub$topic <- jtopic
  return(jsub)
}) %>%
  bind_rows()

# Save to output  ---------------------------------------------------------



# save to output if not already
obj.outf <- "~/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
if (!file.exists(obj.outf)){
  save(out.objs, file = obj.outf)
}

# Plot the 4 umaps with proper colors  ------------------------------------

jsub <- dat.umap.long.trajs[[1]]
jcolor <- colsvec[[1]]

m <- PlotXYNoColor(jsub, xvar = "umap1", yvar = "umap2", jcol = jcolor, jsize = 1)
print(m)

mlst <- lapply(jmarks, function(jmark) PlotXYNoColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", jcol = colsvec[[jmark]], jsize = 1))




# Plot topics -------------------------------------------------------------

# eryth, bcells, nk cells, granu, progenitors

# topics.vec <- c(7, 14, 40, 47, 48)  
topics.vec <- c(1, 44, 24, 50, 11, 16, 19)  
hits <- c("Hbb-bs", "Il2ra", "Prf1", "S100a8", "Kit", "Gzmb", "Inpp4b", "Irf8", "S100a7a", "S100a6", "Tal1", "Sox6", "Gypa", "Car2", "Junb", "Atf3", "Runx2", "Etv6", "Sirpa", "F13a1", "Il10ra", "Lect2")
# hits <- c("Hbb-bs", "Il2ra")
topname <- c("Erythroblasts", "B cells", "NK cells", "Granulocytes", "Progenitors")

 
# Plot hits for H3K4me1 ---------------------------------------------------

topics.mat <- as.data.frame(posterior(out.objs[[jmark]]$out.lda)$topics)
colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))

terms.keep <- subset(terms.filt, gene %in% hits)$term
imput.mat <- t(posterior(out.objs[[jmark]]$out.lda)$topics %*% posterior(out.objs[[jmark]]$out.lda)$terms[, terms.keep])  # not enough memory?!?!?
imput.sub <- tidyr::gather(data.frame(term = terms.keep, imput.mat[terms.keep, ]), key = "cell", value = "exprs", -term)
imput.sub$cell <- gsub("\\.", "-", imput.sub$cell)
imput.sub <- left_join(imput.sub, terms.filt)
imput.wide <- tidyr::spread(imput.sub %>% dplyr::select(gene, exprs, cell), key = "gene", value = "exprs")

topics.mat$cell <- rownames(topics.mat)

# for topics and gene
dat.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat %>% dplyr::select(c("cell", paste0("topic_", topics.vec))))
dat.merge <- left_join(dat.merge, imput.wide)

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/umaps_and_hits/umaps_and_hits.", jmark, ".pdf"), useDingbats = FALSE)
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)

# Show H3K4me1 with louvain and trajectories ------------------------------

traj.jsize <- 2
ctypes <- c("eryth", "granu", "lymphoid", "mega", "nk")
cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")
m.trajs <- ggplot(dat.umap.long.trajs[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank()) +
  scale_color_manual(values=cbPalette2) + 
  xlab("") + ylab("")
for (ctype in ctypes){
  m.trajs <- m.trajs + geom_path(data = trajs[[jmark]][[ctype]], color = "gray25", size = traj.jsize, alpha = 0.5)
}
print(m.trajs)



m1 <- ggplot(zscores.avg.all %>% filter(topic %in% topics.vec), aes(y = topic, x = CellType, size = weight, color = weight)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m1)

for (jtop in topics.vec){
  m1 <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = paste0("topic_", jtop), jsize = 2, jtitle = paste0("topic_", jtop), jcol = scales::muted("darkred"))
  print(m1)
  
  # add gene exprs
  jsub.avg <- subset(zscores.avg.all %>% filter(topic == jtop)) %>%
    arrange(desc(zscore))
  ctypevec <- jsub.avg$CellType
  jsub <- subset(dat.long, Gene_Name %in% top.genes.all.lst[[as.character(jtop)]])
  jsub.avg$CellType <- factor(jsub.avg$CellType, levels = ctypevec)
  jsub$CellType <- factor(jsub$CellType, levels = ctypevec)
  # reorder by avg?
  # sort things
  m.genes <- ggplot(jsub, aes(x = CellType, y = zscore)) + geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) + 
    ggtitle(paste("Topic", jtop))
  m.top <- ggplot(jsub.avg, aes(x = CellType, y = weight)) + geom_bar(stat = "identity") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) + 
    ggtitle(paste("Topic:", jtop))
  hits.sub <- subset(terms.filt, topic == jtop & rnk <= keep.top.genes)
  hits.sub <- OrderDecreasing(hits.sub, jfactor = "termgene", jval = "weight")
  m.hits <- ggplot(hits.sub, aes(x = termgene, y = log10(weight), label = gene)) +
    geom_text_repel() +
    theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Topic:", jtop))
  print(m.genes)
  # print(m.top)
  print(m.hits)
}

for (hit in hits){
  print(hit)
  if (hit %in% colnames(dat.merge)){
    m1 <- PlotXYWithColor(dat.merge, xvar = "umap1", yvar = "umap2", cname = paste0("`", hit, "`"), jsize = 2, jtitle = hit, strip.ticks = TRUE)
    print(m1)
  } else {
    print(paste("Skipping:", hit))
    next
  }
}
dev.off()

# # Atf3 shows good signal, but not Il10ra. Why? 
# subset(terms.filt, gene == "Atf3")
# 
# load(paste0("/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_", jmark, "_bin_TRUE_k_50.genomewide_nofilt.RData"), v=T)  # terms.filt
# terms.all <- terms.filt

# why is Atf3 also in granulocytes? 
dat.merge.alltopics <- left_join(dat.umap.long.trajs[[jmark]], topics.mat %>% dplyr::select(c("cell", paste0("topic_", seq(50)))))
# jgene <- "Sox6"
# dat.merge.alltopics <- left_join(dat.merge.alltopics, )
jgene <- "Atf3"
jgene <- "Etv6"
jgenes <- c("Etv6", "Atf3", "Sirpa")

print(head(subset(terms.all, gene == jgene)))
jtop <- (terms.all %>% filter(gene == jgene) %>% filter(weight == max(weight)))$topic[[1]]
# plot topic
m1 <- PlotXYWithColor(dat.merge.alltopics, xvar = "umap1", yvar = "umap2", cname = paste0("topic_", jtop), jsize = 2, jtitle = paste0("topic_", jtop), jcol = scales::muted("darkred"))
print(m1)

# plot gene

terms.keep <- subset(terms.filt, gene %in% jgenes)$term
imput.mat <- t(posterior(out.objs[[jmark]]$out.lda)$topics %*% posterior(out.objs[[jmark]]$out.lda)$terms[, terms.keep])  # not enough memory?!?!?
imput.sub <- tidyr::gather(data.frame(term = terms.keep, imput.mat[terms.keep, ]), key = "cell", value = "exprs", -term)
imput.sub$cell <- gsub("\\.", "-", imput.sub$cell)
imput.sub <- left_join(imput.sub, terms.filt)
imput.wide <- tidyr::spread(imput.sub %>% dplyr::select(gene, exprs, cell), key = "gene", value = "exprs")

dat.merge.alltopics <- left_join(dat.merge.alltopics, imput.wide)

jgene <- "Etv6"
jgene <- "Atf3"
m1.gene <- PlotXYWithColor(dat.merge.alltopics, xvar = "umap1", yvar = "umap2", cname = jgene, jsize = 2, jtitle = jgene, jcol = scales::muted("darkred"))
print(m1.gene)






# Jake Yeung
# Date of Creation: 2019-09-30
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/2-quick_SVD_exploration.R
# Exploration of SVD


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)

library(irlba)

library(Matrix)
library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(ggrepel)

# Load data ---------------------------------------------------------------


inf <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells/H3K4me1_PZ-StemCells_matsMergedNonenriched_2019-09-29.RData"

load(inf, v=T)

lsi.out <- RunLSI(as.matrix(count.dat$counts), n.components = 50)

# plot out
jsettings <- umap.defaults
jsettings$n_neighbors <- 25
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(lsi.out$u, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

dat.umap.long$is.stem <- sapply(dat.umap.long$cell, function(x) grepl("stem-cell", x))
ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + 
  facet_wrap(~is.stem) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.Robj"
load(inf.lda, v=T)

out.lda <- out.lda[[3]]

umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
dat.umap.long.lda <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

dat.umap.long.lda$is.stem <- sapply(dat.umap.long.lda$cell, function(x) grepl("stem-cell", x))

ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2)) + geom_point() + 
  facet_wrap(~is.stem) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.long.lda <- DoLouvain(topics.mat = posterior(out.lda)$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.lda)

ggplot(dat.umap.long.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  facet_wrap(~is.stem) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# do some variance argument?

inf.proj <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/projections/H3K4me1_matsMergedNonenriched.matsMergedEnriched.RData"
load(inf.proj, v=T)

dat.pred <- predict(umap.out, out.lda.predict$topics)

dat.pred.long <- data.frame(cell = rownames(dat.pred), umap1 = dat.pred[, 1], umap2 = dat.pred[, 2], stringsAsFactors = FALSE) %>%
  mutate(is.stem = TRUE)

dat.umap.pred.merged <- bind_rows(subset(dat.umap.long.lda, select = -louvain), dat.pred.long)

ggplot(dat.umap.pred.merged, aes(x = umap1, y = umap2)) + geom_point() + 
  facet_wrap(~is.stem) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Look at topics  ---------------------------------------------------------

# where are eryths?
inf.outobj <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.OutObjs.RData"
load(inf.outobj, v=T)

terms.long <- data.frame(term = colnames(out.objs$tm.result$terms), as.data.frame(t(out.objs$tm.result$terms)), stringsAsFactors = FALSE) %>%
  gather(key = "topic", value = "weight", -term) %>%
  mutate(topic = gsub("X", "topic_", topic)) %>%
  group_by(topic) %>%
  arrange(desc(weight)) %>%
  mutate(rnk = seq(length(weight))) %>%
  rowwise()

# add gene name
terms.long <- left_join(terms.long, subset(out.objs$regions.annotated, select = c(region_coord, SYMBOL)), by = c("term" = "region_coord"))

topics.sum <- OrderTopicsByEntropy(out.objs$tm.result) %>%
  mutate(topic = gsub("X", "topic_", topic))

# plot topic loadings onto the data
colnames(out.objs$tm.result$topics) <- paste0("topic_", colnames(out.objs$tm.result$topics))
colnames(out.lda.predict$topics) <- paste0("topic_", colnames(out.lda.predict$topics))

topics.mat.merged <- rbind(out.objs$tm.result$topics, out.lda.predict$topics)
topics.mat.merged <- data.frame(cell = rownames(topics.mat.merged), rbind(out.objs$tm.result$topics, out.lda.predict$topics))

# plot output
dat.umap.pred.merged.joined <- left_join(dat.umap.pred.merged, topics.mat.merged)

jtopic <- topics.sum$topic[[1]]
PlotXYWithColor(dat.umap.pred.merged.joined, xvar = "umap1", yvar = "umap2", cname = jtopic) + scale_color_viridis_c() + facet_wrap(~is.stem)

# what are the genes?
print(data.frame(subset(terms.long, topic == "topic_15")))


# Load public data  -------------------------------------------------------


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  group_by(Gene_Name, CellType) %>%
  summarise(FPKM = sum(FPKM)) %>%
  rowwise() %>%
  mutate(logFPKM = log2(FPKM + 1))

# normalize across samples?
ggplot(dat.long, aes(x = CellType, y = logFPKM)) + geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dat.mat <- tidyr::spread(dat.long %>%
                           ungroup() %>%
                           # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>%
                           mutate(gene = Gene_Name) %>%
                           dplyr::select(gene, CellType, logFPKM),
                         key = CellType, value = logFPKM)  %>%
  as.data.frame()
rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL

cnames.tmp <- colnames(dat.mat)
rnames.tmp <- rownames(dat.mat)
dat.mat <- preprocessCore::normalize.quantiles(as.matrix(dat.mat), copy = TRUE)  # strong normalization,
colnames(dat.mat) <- cnames.tmp
rownames(dat.mat) <- rnames.tmp

boxplot(dat.mat)

dat.long.norm <- tidyr::gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "logFPKMnorm", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(logFPKMnorm, center = TRUE, scale = TRUE))


# Get top topics ----------------------------------------------------------

keeptop <- 100
jtopic <- "topic_15"

pdf(paste0("/Users/yeung/data/scchic/pdfs/stemcell_analysis/stemcell_analysis.", Sys.Date(), ".pdf"))

for (jtopic in topics.sum$topic){
  print(jtopic)
  jsub.terms <- subset(terms.long, topic == jtopic & rnk <= keeptop)
  jgenes <- jsub.terms$SYMBOL
  
  m.umap <- PlotXYWithColor(dat.umap.pred.merged.joined, xvar = "umap1", yvar = "umap2", cname = jtopic) + scale_color_viridis_c() + facet_wrap(~is.stem)
  m.genes <- ggplot(subset(dat.long.norm, gene %in% jgenes), aes(x = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE), y = zscore)) + 
    geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggtitle(jtopic)
  # plot top genes
  m.top <- jsub.terms %>%
    ggplot(aes(x = forcats::fct_reorder(.f = term, .x = weight, .fun = median, .desc = TRUE), y = log10(weight), label = SYMBOL)) +
    geom_point(size = 0.25) +
    theme_bw(8) +
    geom_text_repel(size = 3.5, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(paste("Top peak weights for:", jtopic))
  print(m.umap)
  print(m.genes)
  print(m.top)
}

dev.off()


# Do intrachromosomal -----------------------------------------------------

# topics.mat.all <- rbind(out.objs$tm.result$topics, out.lda.predict$topics)
# terms.mat <- out.objs$tm.result$terms

# do variance across chromosomes
jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
jchromos.grep <- paste(jchromos, ":", sep = "")
#jmark <- "H3K4me1"
# jmarks <-

jfac <- 10^6
jpseudo <- 0
dat.impute.log <- log2(t(out.objs$tm.result$topics %*% terms.mat) * jfac + jpseudo)
dat.impute.log.enrich <- log2(t(out.lda.predict$topics %*% out.lda.predict$terms) * jfac + jpseudo)

cells.sd <- GetCellSd(dat.impute.log, "", log2.scale = FALSE, fn = var)

cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)
cells.var.chromo$experi <- "Control"
cells.var.chromo.enrich <- CalculateVarAll(dat.impute.log.enrich, jchromos)
cells.var.chromo.enrich$experi <- "StemCellEnriched"

cells.var.chromo.merge <- rbind(cells.var.chromo, cells.var.chromo.enrich)

umap.out.long.merge <- left_join(dat.umap.pred.merged.joined, cells.var.chromo.merge)

# umap.out.long$cell.var.within.sum.norm <- umap.out.long$cell.var.within.sum.norm
m.intravar <- PlotXYWithColor(umap.out.long.merge, xvar = "umap1", yvar = "umap2", cname = "cell.var.within.sum.norm") + facet_wrap(~experi) + scale_color_viridis_c(direction = -1)

# do enrichment of variance
m.density1 <- ggplot(umap.out.long.merge, aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density() +
  theme_bw(24)+ theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank()) +
  facet_wrap(~experi, ncol = 1) +
  xlab("Intra-chromosomal Variance")
m.density2 <- ggplot(umap.out.long.merge, aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(position = "dodge", alpha = 0.5) +
  theme_bw(24)+ theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank()) +
  xlab("Intra-chromosomal Variance")

pdf(paste0("/Users/yeung/data/scchic/pdfs/stemcell_analysis/stemcell_analysis.variance.", Sys.Date(), ".pdf"))
  print(m.intravar)
  print(m.density1)
  print(m.density2)
dev.off()
  

  





# Jake Yeung
# Date of Creation: 2019-04-21
# File: ~/projects/scchic/scripts/scripts_analysis/tf_activity_debugging/explore_TSS_LDA_with_binned_umap.R
# If we use TSS LDA, how does it look like on the umap? 


rm(list=ls())

library(topicmodels)
library(dplyr)
library(ggplot2)
library(tidytext)
library(umap)
library(data.table)
library(tidyr)

jscale <- 10^7
jpseudo <- 0
.log <- FALSE

jmark <- "H3K4me1"

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")



# Load UMAP ---------------------------------------------------------------

inf.umap <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.umap))
load(inf.umap, v=T)

# Load TSS LDA ------------------------------------------------------------

tssdist <- 80000
jdate <- "2019-04-20"
Kvec <- "50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, "_GeneTSS.Dedup.", jdate, ".", tssdist, ".Robj")
# inf.lda <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-50_GeneTSS.Dedup.2019-04-21.20000.Robj"

assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)


# Load bulk ---------------------------------------------------------------

dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))


# Plot UMAP, overlay genes ------------------------------------------------


out.lda <- out.lda[[length(out.lda)]]

# cnames.old <- unname(colnames(mat.impute))
# cnames.new <- SwitchColnames(cnames.old, jsplit = "-")
# colnames(mat.impute) <- cnames.new

out.lda@documents <- SwitchColnames(unname(out.lda@documents), jsplit = "-")

tm.result <- posterior(out.lda)

# get UMAP coords
umap.out <- umap(tm.result$topics)

umap.long <- data.frame(cell = unname(rownames(umap.out$layout)), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)


top.cells <- tidy(out.lda, matrix = "gamma") %>%
  group_by(topic) %>%
  arrange(desc(gamma)) %>%
  mutate(rnk = seq(length(gamma))) %>%
  mutate(gamma.zscore = scale(gamma, center = TRUE, scale = TRUE)) %>%
  dplyr::rename(cell = document)

top.cells.sum <- top.cells %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(gamma.zscore < quantile(gamma.zscore, 0.95)) %>%
  mutate(zscore.prob = exp(gamma.zscore) / sum(exp(gamma.zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)

top.peaks <- tidytext::tidy(out.lda, matrix = "beta", log = FALSE) %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta))) %>%
  mutate(beta.zscore = scale(beta, center = TRUE, scale = TRUE)) %>%
  rowwise() %>%
  mutate(gene = strsplit(term, ";")[[1]][[2]])


mat.impute <- t(tm.result$topics %*% tm.result$terms)

# filter top terms
m.celldens <- ggplot(top.cells, aes(x = gamma.zscore)) + geom_density() + facet_wrap(~topic) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.celldens)

m.umap <- ggplot(umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m.umap)

m.umap.orig <- ggplot(dat.umap.long.trajs[[jmark]], aes(x = umap1, y = umap2)) + geom_point()

print(m.umap.orig) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Can we plot top genes ---------------------------------------------------


# get gene list
rnames <- rownames(mat.impute)
rnames.keep <- grepl(";", rnames)


mat.impute.sub <- mat.impute[rnames.keep, ]

genes <- sapply(rownames(mat.impute.sub), function(x){
  g <- tryCatch({
    return(strsplit(x, ";")[[1]][[2]])
  }, error = function(e){
    return("Peak")
  })
}, USE.NAMES = FALSE)

exprs.long <- data.frame(peak = rownames(mat.impute.sub), gene = genes, as.data.frame(mat.impute.sub)) %>%
  tidyr::gather(key = "cell", value = "exprs", c(-peak, -gene))

exprs.long <- left_join(exprs.long, dat.umap.long.trajs[[jmark]] %>% dplyr::select(cell, umap1, umap2, louvain))

jgenes <- c("Hbb-bs", "Gata1", "Foxo1", "Inpp4b", "S100a8", "Hs3st5", "Il2ra", "Prf1", "Sox6", "Gata2", "Pax5", "Ly6c2", "Gypa")
jgenes <- c(jgenes, "Tal1", "Mbd2", "Bcl3", "Foxc1", "Nrf1", "Hmbox1", "Spi1", "Gfi1", "Ebf3", "Cebpd", "Cebpb", "Cebpg", "Cebpa", "Pax6", "Pou2f2", "Ebf1")

pdf(paste0("/tmp/genes_H3K4me1_TSS_on_umap.", ".", tssdist, ".", jdate, ".", Kvec, ".pdf"), useDingbats = FALSE)
for (jgene in jgenes){
  print(jgene)
  jsub <- subset(exprs.long, gene == jgene)
  jsub.annot <- subset(top.peaks, gene == jgene)
  if (nrow(jsub) == 0){
    print(paste("Skipping", jgene))
    next
  }
  jpeak <- jsub.annot$term[[1]]
  m1 <- PlotXYWithColor(jsub, xvar = "umap1", yvar = "umap2", cname = "exprs", jtitle = paste(jgene, jpeak))
  print(m1)
}
dev.off()


# Find NK cells? ----------------------------------------------------------
    
    


# plot interestheing topics
top.cells.with.umap <- left_join(top.cells, dat.umap.long.trajs[[jmark]] %>% dplyr::select(cell, umap1, umap2, louvain))

print(m.celldens)
print(top.cells.sum)

jtops <- top.cells.sum$topic

topn <- 200

pdf(paste0("/tmp/topics_H3K4me1_20kb_TSS_on_umap.", topn, ".", tssdist, ".", jdate, ".", Kvec, ".pdf"), useDingbats = FALSE)
  for (jtop in jtops){
    print(jtop)
    m1 <- PlotXYWithColor(top.cells.with.umap %>% filter(topic == jtop), xvar = "umap1", yvar = "umap2", cname = "gamma") + ggtitle(jtop)
    print(m1)
    
    top.genes <- subset(top.peaks, topic == jtop)$gene[1:topn]
    jsub <- subset(dat.long, Gene_Name %in% top.genes)
    jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
    jlevels <- as.character(jsub.sorted.summarised$CellType)
    jsub$CellType <- factor(jsub$CellType, levels = jlevels)
    m.bulk <- ggplot(jsub,
           aes(x = CellType , y = zscore)) +
      geom_boxplot() +
      # geom_violin() +
      geom_jitter(width = 0.1, size = 0.5) +
      # geom_line() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste0("Topic", jtop))
    print(m.bulk)
  }
dev.off()

# list top hits
jtopics <- c(34, 41, 40)

for (jtopic in jtopics){
  print(subset(top.peaks, topic == jtopic))
}

# erythry genes
head(as.data.frame(subset(top.peaks, topic == 17)), n = 150)
head(as.data.frame(subset(top.peaks, topic == 9)), n = 20)

jsub <- subset(dat.long, Gene_Name == "Igkv12-44") %>%
  arrange(desc(zscore))
jsub$CellType <- factor(jsub$CellType, levels = jsub$CellType)
ggplot(jsub, aes(x = CellType, y = zscore)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(hjust = 1,vjust = 1,angle=45))

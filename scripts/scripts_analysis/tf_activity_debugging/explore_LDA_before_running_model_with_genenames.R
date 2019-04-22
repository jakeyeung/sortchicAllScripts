# Jake Yeung
# Date of Creation: 2019-04-20
# File: ~/projects/scchic/scripts/scripts_analysis/tf_activity_debugging/explore_LDA_before_running_model_with_genenames.R
# Explore lDA for terms that have gene names at end

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

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load UMAPs and LDA ------------------------------------------------------

# jsuffix <- "_GeneTSS.Dedup"
jsuffix <- "_GeneTSS.Dedup"
# jsuffix <- "_GeneTSS.Dedup.RbindHiddenDomains"
jdate <- "2019-04-20"
tssdist <- 20000
Kvec <- "25"
inf.spring <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
# inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-25_50.Robj"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, jsuffix, ".", jdate, ".", tssdist, ".Robj")

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.spring))
load(inf.spring, v=T)
load(inf, v=T)

out.lda <- out.lda[[length(out.lda)]]

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
  filter(gamma.zscore < quantile(gamma.zscore, 0.98)) %>%
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


# top.peaks.sum <- top.peaks %>%
#   group_by(term) %>%
#   mutate(topics.zscore = scale(beta, center = TRUE, scale = TRUE)) %>%
#   mutate(topics.zscore.prob = exp(topics.zscore) / sum(exp(topics.zscore))) %>%
#   summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
#   arrange(entropy)

mat.impute <- t(tm.result$topics %*% tm.result$terms)

# filter top terms
m.celldens <- ggplot(top.cells, aes(x = gamma.zscore)) + geom_density() + facet_wrap(~topic) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.celldens)

m.umap <- ggplot(umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(m.umap)

# get gene list
rnames <- rownames(mat.impute)
rnames.keep <- grepl(";", rnames)

# find tal1 peak
# rnames.tal1 <- grepl("chr4:1150[5-6][0-9]{4}", rnames)
# rnames.keep <- rnames.keep | rnames.tal1

cnames.old <- unname(colnames(mat.impute))
cnames.new <- SwitchColnames(cnames.old, jsplit = "-")
colnames(mat.impute) <- cnames.new
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

# Plot UMAP and overlay gene expression -----------------------------------

exprs.long <- left_join(exprs.long, dat.trajs.long)

# # jgene <- "Ebf1"
# jgene <- "T"
# jsub <- subset(exprs.long, peak == "chr4:115065000-115066000") %>%
#   mutate(exprs = log10(exprs * jscale + jpseudo))
# jpeak <- jsub$peak[[1]]
# m1 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "exprs", jtitle = paste(jgene, jpeak))
# print(m1)

# 
jgenes <- c("Hbb-bs", "Gata1", "Foxo1", "Inpp4b", "S100a8", "Hs3st5", "Il2ra", "Prf1", "Klf", "Sox6", "Gata2", "Pax5", "Ly6c2", "Gypa")
jgenes <- c(jgenes, "Tal1", "Mbd2", "Bcl3", "Foxc1", "Nrf1", "Hmbox1", "Spi1", "Gfi1", "Ebf3", "Cebpd", "Cebpb", "Cebpg", "Cebpa", "Pax6", "Pou2f2", "Ebf1")

pdf(paste0("/tmp/lda_check.log.", .log, ".", jsuffix, ".", tssdist,  ".", jdate, ".pdf"), useDingbats = FALSE)
for (jgene in jgenes){
  if (.log){
    jsub <- subset(exprs.long, gene == jgene) %>%
      mutate(exprs = log2(exprs * jscale + jpseudo))
  } else {
    jsub <- subset(exprs.long, gene  == jgene)
  }
  if (nrow(jsub) == 0){
    print(paste("Skipping", jgene))
    next
  }
  jpeak <- jsub$peak[[1]]
  m1 <- PlotXYWithColor(jsub, xvar = "X1", yvar = "X2", cname = "exprs", jtitle = paste(jgene, jpeak))
  print(m1)
}
dev.off()



# What are the top hits? --------------------------------------------------

# which peaks are most interesting?


# plot interestheing topics
top.cells <- left_join(top.cells, umap.long)

print(m.celldens)

jtop <- 1  # bcells
jtop <- 4  # granulocytes
jtop <- 21  # eryths

jtop <- 43
jtop <- 28

jtop <- 44
jtop <- 8
jtop <- 46

jtop <- 39
jtop <- 48
jtop <- 9
jtop <- 37
jtop <- 15
m1 <- PlotXYWithColor(top.cells %>% filter(topic == jtop), xvar = "umap1", yvar = "umap2", cname = "gamma") + ggtitle(jtop)
print(m1)


# Do top genes look OK comparing with bulk?  ------------------------------

dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

jtopic <- 28
jtopic <- 21
jtopic <- 43

jtopic <- 44
jtopic <- 8
jtopic <- 46

jtopic <- 39
jtopic <- 48
jtopic <- 15

jtopic <- 50
print(data.frame(head(subset(top.peaks, topic == jtopic), n = 20)))

topn <- 100
top.genes <- subset(top.peaks, topic == jtopic)$gene[1:topn]
jsub <- subset(dat.long, Gene_Name %in% top.genes)
jsub.sorted.summarised <- jsub %>% group_by(CellType) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(CellType)
jlevels <- as.character(jsub.sorted.summarised$CellType)
jsub$CellType <- factor(jsub$CellType, levels = jlevels)
ggplot(jsub,
       aes(x = CellType , y = zscore)) +
  geom_boxplot() +
  # geom_violin() +
  geom_jitter(width = 0.1, size = 0.5) +
  # geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(jtopic)


top.genes[which(top.genes %in% jsub$Gene_Name)]

jgene <- "Gm21738"
jgene <- "Lalba"
jgene <- "Gm10801"
jgene <- "Gm10800"
jgene <- "Gm27940"
jgene <- "Gm44385"
jgene <- "Igkv1-110"
jgene <- "Kat2b"
jgene <- "Dennd2c"
jgene <- "Tnpo1"
jgene <- "Pcx"
jgene <- "Tal1"

jgene <- "Cebpb"
ggplot(dat.long %>% filter(Gene_Name == jgene), aes(x = CellType, y = zscore)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(jgene)

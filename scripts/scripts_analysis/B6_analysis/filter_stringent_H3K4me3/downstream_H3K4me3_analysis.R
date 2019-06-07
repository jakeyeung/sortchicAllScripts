# Jake Yeung
# Date of Creation: 2019-06-04
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/H3K4me3_stringent_analysis/downstream_H3K4me3_analysis.R
# Downstream analysis of H3K4me3 with more stringent cutoffs 
# 

rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(umap)
library(igraph)
library(hash)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)


source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- "H3K4me3"
jbin <- TRUE
keep.top.genes <- 150

# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))


# Load LDA ----------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.stringent_filter/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50.Robj"
assertthat::assert_that(file.exists(inf))

kchoose <- 50
out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)

# out.lda <- out.objs$out.lda

# load(inf, v=T)

# # take K=50
# out.lda <- out.lda[[4]]
# assertthat::assert_that(out.lda@k == 50)


# Get posteriors ----------------------------------------------------------

tm.result <- posterior(out.objs$out.lda)

# Do UMAP  ----------------------------------------------------------------

mindist <- 0.4
nneigh <- 34
randomstate <- 123
jmetric <- "euclidean"

# umap.settings <- GetUmapSettings(nn = nneigh, jmindist = mindist, jmetric = jmetric, seed = randomstate)
umap.settings <- umap.defaults
umap.settings$min_dist <- mindist; umap.settings$n_neighbors <- nneigh

umap.out <- umap(tm.result$topics, config = umap.settings)

dat.umap.long <- data.table(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  mutate(umap2 = umap2 * -1)  # flip y axis 

PlotXYNoColor(dat.umap.long, xvar = "umap1", yvar = "umap2")


# Check FACS granulocytes -------------------------------------------------

dat.facs.filt <- lapply(jmarks, LoadFACSGetLoadings) %>%
  bind_rows()

dat.facs.merge <- left_join(dat.umap.long, dat.facs.filt %>% filter(mark == "H3K4me3"))

PlotXYWithColor(dat.facs.merge, xvar = "umap1", yvar = "umap2", cname = "loadings", jsize = 5)

# Annotated and find relevant genes  ---------------------------------------


# Do louvain clustering ---------------------------------------------------

nn.louv <- 34
jmetric.louv <- "euclidean"
jmindist.louv <- 0.4
jsettings <- GetUmapSettings(nn.louv, jmetric.louv, jmindist.louv, randomstate)

louv.hash <- DoLouvain(tm.result$topics, jsettings)

# add to umap
dat.umap.long$louvain <- as.character(sapply(as.character(dat.umap.long$cell), function(x) louv.hash[[x]]))

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
# cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")

m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette) + ggtitle(paste(nn.louv))


# Find genes --------------------------------------------------------------

topics.mat <- data.frame(cell = rownames(tm.result$topics), as.data.frame(tm.result$topics))
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

terms.long <- data.frame(term = colnames(tm.result$terms), as.data.frame(t(tm.result$terms)), stringsAsFactors = FALSE) %>%
  gather(key = "topic", value = "weight", -term) %>%
  mutate(topic = gsub("X", "", topic)) %>%
  group_by(topic) %>%
  arrange(desc(weight)) %>%
  mutate(rnk = seq(length(weight))) %>%
  rowwise()

# Load facs ---------------------------------------------------------------



# Reannotate genes --------------------------------------------------------


# filter out top 1000 peaks

# annotate terms
terms.filt.top <- terms.long %>%
  # filter(rnk < 1000) %>%  # DO GENOME WIDE
  rowwise()

tss.dat <- fread("/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", col.names = c("seqnames", "start", "end", "tssname"))
tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, ";")[[1]][[2]])

annots.biomart <- out.objs$regions.annotated %>% 
  mutate(midpt = start + (end - start) / 2) %>%
  filter(region_coord %in% terms.filt.top$term)


annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% dplyr::select(seqnames, start, end, SYMBOL, region_coord), keep.extra.columns = TRUE)

annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)

out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")

out2.df = data.frame(annots.gr[queryHits(out2),], annots.tss.gr[subjectHits(out2),]) %>%
  mutate(midpt = start + round(width / 2),
         midpt.1 = start.1 + round(width.1 / 2),
         dist.to.tss = midpt.1 - midpt)

# filter closest
out2.df.closest <- out2.df %>%
  group_by(region_coord) %>%
  filter(abs(dist.to.tss) == min(abs(dist.to.tss)))

# Find interesting topics -------------------------------------------------

terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, sep = ";")
terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)  # error stack overflow?

terms.filt <- terms.filt.top %>%
  mutate(termgene = ifelse(!is.null(terms.hash[[term]]), terms.hash[[term]], NA)) %>%
  filter(!is.na(termgene)) %>%
  mutate(gene = sapply(termgene, function(x) strsplit(x, ";")[[1]][[2]])) %>%
  group_by(gene)  # dont filter use 
# filter(weight == max(weight))




# plot --------------------------------------------------------------------


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



# Order topics by celltype identity ---------------------------------------

zscores.sum <- zscores.avg.all %>%
  group_by(topic) %>%
  summarise(entropy = -sum(weight * log(weight))) %>%
  arrange(entropy)


# Plot top hits -----------------------------------------------------------

print(topics.sum)
# plot a topic
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X24")
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X13")
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X26")
# PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = "X34")

ntopicsfilt <- round(nrow(topics.sum) / 2)
# ntopicsfilt <- nrow(topics.sum)
# topics.keep <- topics.sum$topic[1:ntopicsfilt]
topics.keep <- zscores.sum$topic[1:ntopicsfilt]
names(topics.keep) <- topics.keep

topics.ordered <- zscores.sum$topic

genes.keep <- subset(terms.filt, topic %in% topics.keep)$gene # remove meaningless topics



top.genes.lst <- lapply(topics.keep, function(jtopic){
  jtmp <- subset(terms.filt, topic == jtopic)$gene
  jtmp <- jtmp[1:min(length(jtmp), keep.top.genes)]
  return(jtmp)
})

# do all topics



# get averag zscore across celltypes for each topic
zscores.avg <- lapply(topics.keep, function(jtopic){
  jtopic <- as.character(jtopic)
  top.genes <- top.genes.lst[[jtopic]]
  jsub <- subset(dat.long, Gene_Name %in% top.genes) %>%
    group_by(CellType) %>%
    summarise(zscore = mean(zscore)) %>%
    ungroup() %>%
    mutate(weight = exp(zscore) / sum(exp(zscore)))
  jsub$topic <- jtopic
  return(jsub)
}) %>%
  bind_rows()


# Check top topics --------------------------------------------------------

# topvecs <- topics.sum$topic[1:10]

dat.withtopic <- left_join(dat.umap.long, topics.long) 

pdf(paste0("~/data/scchic/pdfs/B6_figures/B6_", jmark, "_bin_", jbin, "_k_", kchoose, "_topics_plot_", Sys.Date(), ".stringent_filter.pdf"), useDingbats = FALSE)
for (topvec in topics.keep){
  m1 <- PlotXYWithColor(dat.withtopic %>% filter(topic == topvec), xvar = "umap1", yvar = "umap2", cname = "weight", jtitle = topvec)
  print(m1)
}
dev.off()




# do all the topics 

# this goes into main figure probably

pdf(paste0("~/data/scchic/pdfs/B6_figures/B6_", jmark, "_bin_", jbin, "_k_", kchoose, "_celltypes_by_topic_only.keepgenes.", keep.top.genes, ".", Sys.Date(), ".winbugfix.stringent_filter.pdf"), useDingbats = FALSE)
# ctypes.filt <- c("T_cell", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "granulocyte", "Kit_and_Sca1-positive_hematopoietic_stem_cell")
topskeep <- topics.sum$topic[1:5]
m1 <- ggplot(zscores.avg %>% filter(topic %in% topskeep), aes(y = topic, x = CellType, size = weight, color = weight)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m1)
for (jtop in topics.keep){
  # plot also the UMAP 
  m.umap <- PlotXYWithColor(left_join(dat.umap.long, topics.mat), xvar = "umap1", yvar = "umap2", cname = paste0("X", jtop), jtitle = jtop)
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
  # show top loadings
  hits.sub <- subset(terms.filt, topic == jtop & rnk <= keep.top.genes)
  hits.sub <- OrderDecreasing(hits.sub, jfactor = "termgene", jval = "weight")
  m.hits <- ggplot(hits.sub, aes(x = termgene, y = log10(weight), label = gene)) +
    geom_text_repel() +
    theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Topic:", jtop))
  print(m.umap)
  print(m.genes)
  print(m.top)
  print(m.hits)
}
dev.off()

# Write objects -----------------------------------------------------------

save(terms.filt, out2.df.closest, file = paste0("~/data/scchic/robjs/B6_objs/terms_filt_", jmark, "_bin_", jbin, "_k_", kchoose, ".genomewide_nofilt.stringent_filter.withDist.RData"))

# print(m)

# Save objects and write files  -------------------------------------------

outdir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent"
dir.create(outdir)

pdf(file.path(outdir, paste0("umap_cluster_IDs.", jmark, ".stringent_filter.pdf")))
print(m.louv)
dev.off()

dat.sub <- subset(dat.umap.long, select = c(cell, louvain))
dat.sub$fname <- sapply(dat.sub$cell, function(x) paste0(x, ".sorted.bam"))
for (jclst in unique(base::sort(dat.sub$louvain))){
  dat.subsub <- subset(dat.sub, louvain == jclst)
  outf.sub <- file.path(outdir, paste0(paste("bamlist", jclst, jmark, sep = "-"), ".txt"))
  if (file.exists(outf.sub)){
    next
  } else {
    data.table::fwrite(dat.subsub %>% dplyr::select(fname), file = outf.sub, sep = "\t", col.names = FALSE)
  }
}

save(out.objs, dat.umap.long, umap.settings, jsettings, file = file.path(outdir, paste0("dat_umap_long_with_louvain.", jmark, ".stringent_filter.RData")))

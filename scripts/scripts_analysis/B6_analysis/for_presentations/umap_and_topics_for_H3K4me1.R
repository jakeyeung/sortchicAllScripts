# Jake Yeung
# Date of Creation: 2019-05-18
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/for_presentations/umap_and_topics_for_H3K4me1.R
# Plot pretty for presentation

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(ggrepel)
library(JFuncs)

source("scripts/Rfunctions/PlotFunctions.R")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")

# Load bulk data ----------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
assertthat::assert_that(file.exists(inf.bulkdat))
dat <- fread(inf.bulkdat, sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)


# Load LDA ----------------------------------------------------------------

inf.lda <- "/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
load(inf.lda, v=T)

inf.terms <- "/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_H3K4me1_bin_TRUE_k_50.RData"
load(inf.terms, v=T)  # H3k4me1 only 


pdf(file = paste0("/Users/yeung/data/scchic/pdfs/B6_figures/for_presentation/umap_and_topics.", Sys.Date(), ".pdf"), useDingbats = FALSE)

# Plot stuff --------------------------------------------------------------

# plot H3K9me3 for Fig 4

jmark <- "H3K4me1"
PlotXYNoColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", jcol = "gray80", jsize = 3)

# show topic loadings
topics.mat <- data.frame(cell = rownames(out.objs[[jmark]]$tm.result$topics), out.objs[[jmark]]$tm.result$topics, stringsAsFactors = FALSE)

dat.merge <- left_join(dat.umap.long.trajs[[jmark]], topics.mat)

# replace X with Topic_ in colnames
colnames(dat.merge) <- gsub("^X", "Topic_", colnames(dat.merge))


PlotXYWithColor(dat.merge %>% mutate(Topic_7_norm = scale(Topic_7, center = FALSE, scale = FALSE)), 
                xvar = "umap1", yvar = "umap2", cname = "Topic_7_norm", jsize = 5, 
                jcol = scales::muted('darkred'), jcol.mid = "lightblue", jcol.low = scales::muted('darkblue'), manual.mid = 0)


# Show topic loadings ------------------------------------------------------

# for Erythryoblasts 
jtop <- 7
loadings <- subset(terms.filt, topic == jtop)
loadings <- OrderDecreasing(loadings, "term", "weight")

ggplot(loadings %>% filter(rnk < 70), aes(x = term, y = log10(weight), label = gene)) + 
  theme_bw() + 
  theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_point() + geom_text_repel() + xlab("") + ylab("log10(Region Loading)")

genes.show <- c("Cpeb3", "Sox6", "Tal1", "Aqp", "Hbq", "Hbb", "Slc", "Lmo", "Gyp")
genes.show.grep <- paste(genes.show, collapse = "|")
ggplot(loadings %>% filter(rnk < 70) %>% mutate(gene.lab = ifelse(grepl(genes.show.grep, gene), gene, NA)), aes(x = term, y = log10(weight), label = gene.lab)) + 
  theme_bw() + 
  theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_point() + geom_text_repel(size = 7) + xlab("") + ylab("log10(Region Loading)") + ggtitle("Top 70 regions in erythroblast topic")

ggplot(loadings %>% filter(rnk < 70) %>% mutate(gene.lab = ifelse(grepl(genes.show.grep, gene), gene, NA)), aes(x = term, y = log10(weight), label = gene)) + 
  theme_bw() + 
  theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_point() + geom_text_repel() + xlab("") + ylab("log10(Region Loading)") + ggtitle("Top 70 regions in erythroblast topic")



# Show the bulk data  -----------------------------------------------------

jgenes.keep <- (loadings %>% filter(rnk < 70))$gene

dat.sub <- dat.long %>% filter(Gene_Name %in% jgenes.keep)

# shorten names
dat.sub$CellType <- gsub(pattern = "nucleate_erythrocyte", "erythroblasts", dat.sub$CellType)
dat.sub$CellType <- gsub(pattern = "megakaryocyte-erythroid_progenitor_cell", "erythroid progenitor", dat.sub$CellType)
dat.sub$CellType <- gsub(pattern = "Kit_and_Sca1-positive_hematopoietic_stem_cell", "kit_sca1_pos", dat.sub$CellType)

dat.sub.sum <- dat.sub %>%
  group_by(CellType) %>%
  summarise(zscore = median(zscore)) %>%
  arrange(desc(zscore))
jlevs <- dat.sub.sum$CellType
dat.sub$CellType <- factor(as.character(dat.sub$CellType), levels = jlevs)
ggplot(dat.sub, aes(x = CellType, y = zscore)) + 
  geom_boxplot(outlier.size = 0) + geom_point(size = 0.5) +
  theme_bw() +  
  theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("Cell Type") + ylab("Bulk RNA-seq Levels (Zscore)") + 
  ggtitle("Genes associated with top 70 regions")

# ctypes <- c("granu", "eryth", "mega", "lymphoid")
ctypes <- c("granu")
traj.jsize <- 3
# plot louvains

for (jjmark in jmarks){
  m.louv <- ggplot(dat.umap.long.trajs[[jjmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank()) +
    scale_color_manual(values=cbPalette) +
    xlab("") + ylab("")
  for (ctype in ctypes){
    m.louv <- m.louv + geom_path(data = trajs[[jjmark]][[ctype]], color = "black", size = traj.jsize, alpha = 0.5, arrow = arrow(ends = "last"))
  }
  print(m.louv)
}
m.louv.h3k9me3 <- ggplot(dat.umap.long.trajs[["H3K9me3"]], aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 3) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank()) +
  scale_color_manual(values=cbPalette) +
  xlab("") + ylab("")
for (ctype in ctypes){
  m.louv.h3k9me3 <- m.louv.h3k9me3 + geom_path(data = trajs[["H3K9me3"]][[ctype]], color = "black", size = traj.jsize, alpha = 0.5, arrow = arrow(ends = "first"))
}
print(m.louv.h3k9me3)
dev.off()


# 
# # for Granulocytes 
# jtop <- 47
# loadings <- subset(terms.filt, topic == jtop)
# loadings <- OrderDecreasing(loadings, "term", "weight")
# 
# ggplot(loadings %>% filter(rnk < 70), aes(x = term, y = log10(weight), label = gene)) + 
#   theme_bw() + 
#   theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
#   geom_point() + geom_text_repel() + xlab("") + ylab("log10(Region Loading)")
# 
# # genes.show <- c("Cpeb3", "Sox6", "Tal1", "Aqp", "Hbq", "Hbb", "Slc", "Lmo")
# # genes.show.grep <- paste(genes.show, collapse = "|")
# ggplot(loadings %>% filter(rnk < 70), aes(x = term, y = log10(weight), label = gene)) + 
#   theme_bw() + 
#   theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
#   geom_point() + geom_text_repel() + xlab("") + ylab("log10(Region Loading)")
# 


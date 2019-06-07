# Jake Yeung
# Date of Creation: 2019-06-06
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/transcriptome_versus_marks_topics.R
# Compare transcriptomes: topics analysis


rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)
library(ggrepel)

source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load bulk  --------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
assertthat::assert_that(file.exists(inf.bulkdat))
dat <- fread(inf.bulkdat, sep = "\t")
colnames(dat) <- gsub(" ", "_", colnames(dat))
dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))


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

topics.vec <- c(7, 14, 40, 47, 48, 13)

# Load umaps  -------------------------------------------------------------

load(inf.traj, v=T)

load(inf.termsfilt, v=T)  # terms.filt

load(inf.objs, v=T)



topics.dat <- tidyr::gather(data.frame(cell = rownames(posterior(out.objs$H3K4me1$out.lda)$topics), posterior(out.objs$H3K4me1$out.lda)$topics), key = topic, value = topic.weight, -cell) %>%
  mutate(topic = gsub("X", "", topic))


dat.merged <- left_join(dat.umap.long.trajs[[jmark]], topics.dat)

# Plot topic UMAP then its gene enrichment  -------------------------------

ngenes.keep <- 100
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3")

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/celltypes_top_genes_stringent/topics_vs_bulk_rnaseq.", jmark, ".pdf"), useDingbats = FALSE)
for (jtop in topics.vec){
  m1 <- PlotXYWithColor(dat.merged %>% filter(topic == jtop), xvar = "umap1", yvar = "umap2", cname = "topic.weight", jsize = 2.5, jtitle = paste("Topic:", jtop))
  print(m1)
  # get top genes 
  
  loadings <- subset(terms.filt, topic == jtop)
  loadings <- OrderDecreasing(loadings, "term", "weight")
  
  m1 <- ggplot(loadings %>% filter(rnk < ngenes.keep), aes(x = term, y = log10(weight), label = gene)) + 
    theme_bw() + 
    theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    geom_point() + geom_text_repel() + xlab("") + ylab("log10(Region Loading)") + 
    ggtitle(paste("Topic:", jtop))
  print(m1)
  
  jgenes.keep <- (loadings %>% filter(rnk < ngenes.keep))$gene 
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
  m2 <- ggplot(dat.sub, aes(x = CellType, y = zscore)) + 
    geom_boxplot(outlier.size = 0) + geom_point(size = 0.5) +
    theme_bw() +  
    theme(aspect.ratio=0.3, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    xlab("Cell Type") + ylab("Bulk RNA-seq Levels (Zscore)") + 
    ggtitle(paste0("Genes associated with top ", ngenes.keep, " regions"))
  print(m2)
}
dev.off()
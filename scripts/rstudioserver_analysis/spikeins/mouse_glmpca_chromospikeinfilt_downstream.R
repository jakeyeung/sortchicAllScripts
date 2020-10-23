# Jake Yeung
# Date of Creation: 2020-08-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse_LDA_projections_merged_downstream.R
# 





rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")



hubprefix <- "/home/jyeung/hub_oudenaarden"


jsuffix <- ".chromo2spikeinfilt"

indir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.merged_old_and_new", jsuffix)
dir.create(indir)
inf.meta <- file.path(indir, paste0("cell_cluster_merged_with_spikein_info", jsuffix, ".txt"))

inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein", jsuffix, "/H3K4me3_BM.glmpcaout.penalty_5.RData"))
# inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_mouse_spikein", jsuffix, "/H3K4me3_BM.dev_filt.glmpcaout.penalty_5.RData"))

assertthat::assert_that(file.exists(inf.glmpca))
assertthat::assert_that(file.exists(inf.meta))

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_downstream_BM_H3K4me3/glmpca_downstream.pdf"

pdf(file = outpdf, useDingbats = FALSE)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load orignal data matrix ------------------------------------------------

inf.orig <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda.chromo2spikeinfilt/H3K4me3_BM.match_rownames_with_old.rds"

mat.orig <- readRDS(inf.orig)


# Load meta ---------------------------------------------------------------


dat.meta <- fread(inf.meta)

# Check GLMPCA ------------------------------------------------------------

load(inf.glmpca, v=T)

head(glmpcaout$factors)

dat.umap.glmpca <- DoUmapAndLouvain(glmpcaout$factors, jsettings); dat.umap.glmpca$louvain <- NULL

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

head(spikeincounts)

# inf.spikeins <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_count_tables_mouse_for_lda/H3K4me3_BM.spikeins.txt")
# dat.spikeins <- fread(inf.spikeins)

dat.umap.glmpca <- dat.umap.glmpca %>%
  # left_join(., dat.spikeins, by = c("cell" = "samp")) %>%
  left_join(., subset(dat.meta, select = -c(umap1, umap2)))

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

ggplot(dat.umap.glmpca, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = cbPalette)

dat.dim.glmpca <- data.frame(cell = rownames(glmpcaout$factors), 
                             dim1 = glmpcaout$factors$dim1, 
                             dim2 = glmpcaout$factors$dim2, 
                             dim3 = glmpcaout$factors$dim3, 
                             dim4 = glmpcaout$factors$dim4, 
                             stringsAsFactors = FALSE) %>%
  # left_join(., dat.spikeins, by = c("cell" = "samp")) %>%
  left_join(., dat.meta)

ggplot(dat.dim.glmpca, aes(x = dim1, y = dim2, color = log2(totalcounts / spikeincounts))) +
  geom_point()  + 
  scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.dim.glmpca, aes(x = dim1, y = dim2, color = cluster)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.dim.glmpca, aes(x = dim2, y = dim3, color = cluster)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.meta, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cond) + scale_color_viridis_c(direction = -1)

# plot differences in chromo/spikeincounts 

ggplot(dat.meta %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cond) + scale_color_viridis_c(direction = -1)

ggplot(dat.meta %>% filter(!is.na(spikeincounts) & !cluster %in% c("Dendritic_topic18", "Eryth-Gfi1b_topic7")), aes(x = cluster, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  # facet_wrap(~cond) + 
  scale_color_viridis_c(direction = -1)


# Annotate genes ----------------------------------------------------------

bins <- rownames(glmpcaout$loadings)

inf.tsspretty <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.pretty.bed"
inf.tssorig <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.bed"

dat.tssorig <- fread(inf.tssorig)
dat.tsspretty <- fread(inf.tsspretty)

assertthat::assert_that(identical(dat.tssorig$V1, dat.tsspretty$V1))

tssnamehash <- hash::hash(dat.tssorig$V4, dat.tsspretty$V4)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

bins.annots <- AnnotateCoordsFromList.GeneWise(coords.vec = bins, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

# take top
jdat <- data.frame(bin = rownames(glmpcaout$loadings), glmpcaout$loadings, stringsAsFactors = FALSE) %>%
  arrange(desc(dim2)) %>%
  # arrange(dim2) %>%
  left_join(., bins.annots$out2.df, by = c("bin" = "region_coord"))


# Plot hits on UMAP  ------------------------------------------------------


# inf.tss <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.ByChromo.WithSpikeIns.NoChromo.csv")
inf.mattss <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/mouse/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters/mouse-BM-EtOH-H3K4me3-1.sorted.tagged.countTable.TSS.good_cells_filt.rds")
mat.tss <- readRDS(inf.mattss)
rownames(mat.tss) <- sapply(rownames(mat.tss), function(x) tssnamehash[[x]])


# Plot neutrophil gene ----------------------------------------------------

dat.umap.glmpca.clean <- subset(dat.umap.glmpca, select = c(cell, umap1, umap2))  # for joining
dat.dim.glmpca <- subset(dat.dim.glmpca, select = -c(umap1, umap2))  # for joining
print(head(jdat))
print(tail(jdat))

head(jdat %>% arrange(desc(dim1)), n = 100)

head(bins.annots$out2.df)

tssnames.neutro <- (jdat %>% arrange(desc(dim2)))$tssname[1:10]; tssnames.neutro <- tssnames.neutro[!is.na(tssnames.neutro)]
tssnames.hspcs <- (jdat %>% arrange(dim2))$tssname[1:10]; tssnames.hspcs <- tssnames.hspcs[!is.na(tssnames.hspcs)]
# tssname <- jdat$tssname[[5]]
# subset(jdat, grepl("Hlf", tssname))

# tssname <- "chr3:27292076-27342076;Tnfsf10"
# tssname <- "chr11:90365917-90415917;Hlf"
# tssname <- "chr3:90629301-90679301;S100a7a"
# tssname <- "chr12:118821328-118871328;Sp8"

for (tssname in tssnames.neutro){
  if (is.na(tssname)){
    print("not found, skipping")
    next
  }
  jsub <- data.frame(cuts = mat.tss[tssname, ], cell = colnames(mat.tss), stringsAsFactors = FALSE) %>%
    left_join(., dat.dim.glmpca)  %>%
    left_join(., dat.umap.glmpca.clean)
  
  m1 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(tssname, "Neutro-specific")
  m2 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(tssname, "Neutro-specific")
  print(m1)
  print(m2)
}

for (tssname in tssnames.hspcs){
  if (is.na(tssname)){
    print("not found, skipping")
    next
  }
  jsub <- data.frame(cuts = mat.tss[tssname, ], cell = colnames(mat.tss), stringsAsFactors = FALSE) %>%
    left_join(., dat.dim.glmpca)  %>%
    left_join(., dat.umap.glmpca.clean)
  
  m1 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(tssname, "HSPC-specific")
  m2 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(tssname, "HSPC-specific")
  print(m1)
  print(m2)
}



jsub <- data.frame(cuts = colSums(mat.tss[tssnames.neutro, ]), cell = colnames(mat.tss), stringsAsFactors = FALSE) %>%
  left_join(., dat.dim.glmpca)  %>%
  left_join(., dat.umap.glmpca.clean)

m1 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("Neutros-specific top 10 TSSs")
m2 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("Neutros-specific top 10 TSSs")
print(m1)
print(m2)



jsub <- data.frame(cuts = colSums(mat.tss[tssnames.hspcs, ]), cell = colnames(mat.tss), stringsAsFactors = FALSE) %>%
  left_join(., dat.dim.glmpca)  %>%
  left_join(., dat.umap.glmpca.clean)

m1 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("HSPC-specific top 10 TSSs")
m2 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("HSPC-specific top 10 TSSs")
print(m1)
print(m2)



# jbin <- "chr9:97450000-97500000"

jdat.dim1.high <- jdat %>% arrange(desc(dim1))
jdat.dim1.low <- jdat %>% arrange(dim1)

jbins.dim1.high <- jdat.dim1.high$bin[1:10]
jbins.dim1.low <- jdat.dim1.low$bin[1:10]


for (jbin in jbins.dim1.high){
  
  jsub.orig <- data.frame(cuts = mat.orig[jbin, ], cell = colnames(mat.orig), stringsAsFactors = FALSE) %>%
    left_join(., dat.dim.glmpca) %>%
    left_join(., dat.umap.glmpca.clean)
  
  m1 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(jbin, "dim1 high")
  
  m2 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(jbin, "dim1 high")
  
  print(m1)
  print(m2)
}

jsub.orig <- data.frame(cuts = colSums(mat.orig[jbins.dim1.high, ]), cell = colnames(mat.orig), stringsAsFactors = FALSE) %>%
  left_join(., dat.dim.glmpca) %>%
  left_join(., dat.umap.glmpca.clean)

m1 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("dim1 high sum top 10 bins")

m2 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("dim1 high sum top 10 bins")

print(m1)
print(m2)


for (jbin in jbins.dim1.low){
  
  jsub.orig <- data.frame(cuts = mat.orig[jbin, ], cell = colnames(mat.orig), stringsAsFactors = FALSE) %>%
    left_join(., dat.dim.glmpca) %>%
    left_join(., dat.umap.glmpca.clean)
  
  m1 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(jbin, "dim1 low")
  
  m2 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(jbin, "dim1 low")
  
  print(m1)
  print(m2)
}

jsub.orig <- data.frame(cuts = colSums(mat.orig[jbins.dim1.low, ]), cell = colnames(mat.orig), stringsAsFactors = FALSE) %>%
  left_join(., dat.dim.glmpca) %>%
  left_join(., dat.umap.glmpca.clean)

m1 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("dim1 low sum top 10 bins")

m2 <- ggplot(jsub.orig %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("dim1 low sum top 10 bins")

print(m1)
print(m2)


# Find Grun genes ---------------------------------------------------------

inf.grun <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/bone_marrow_grun/pseudobulk_downsampled_neutrophil_trajectory.2020-05-22.WithTtest.RData"
load(inf.grun, v=T)

jgenes.full <- as.character(ttests.neutroprogs.out$gene[1:20])
jgenes <- sapply(jgenes.full, function(g) strsplit(g, "__")[[1]][[1]])

for (jgene in jgenes){
  jsub <- subset(jdat, gene == jgene)
  if (nrow(jsub) == 0){
    print(paste(jgene, "not found, skipping"))
    next
  }
  tssname <- jsub$tssname[[1]]
  
  jsub <- data.frame(cuts = mat.tss[tssname, ], cell = colnames(mat.tss), stringsAsFactors = FALSE) %>%
    left_join(., dat.dim.glmpca) %>% 
    left_join(., dat.umap.glmpca.clean)
  
  m1 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(tssname, "Grun neut progenitors")
  
  m2 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
    geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c() + ggtitle(tssname, "Grun neut progenitors")
  print(m1)
  print(m2)
}


# add up neutr progenitors
tssnames.all <- subset(jdat, gene %in% jgenes)$tssname

jsub <- data.frame(cuts = colSums(mat.tss[tssnames.all, ]), cell = colnames(mat.tss), stringsAsFactors = FALSE) %>%
  left_join(., dat.dim.glmpca) %>% 
  left_join(., dat.umap.glmpca.clean)

m1 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = dim1, y = dim2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("Grun neut progenitors sum top 20 TSSs")

m2 <- ggplot(jsub %>% filter(!is.na(spikeincounts)), aes(x = umap1, y = umap2, color = log2(cuts / spikeincounts))) + 
  geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle("Grun neut progenitors top 20 TSSs")
print(m1)
print(m2)

dev.off()

# Jake Yeung
# Date of Creation: 2019-11-15
# File: ~/projects/scchic/scripts/scripts_analysis/debug_marker_genes/check_reads_at_TSS_on_UMAP.R
# Check UMAP on TSS

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(GGally)
library(colorRamps)
library(parallel)

library(scchicFuncs)

library(DESeq2)
library(preprocessCore)

library(JFuncs)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
assertthat::assert_that(dir.exists(indir))

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"

winsize <- 50000
indir.raw <- paste0("/Users/yeung/data/scchic/from_cluster/count_mats.fromGeneTSS.0_build95_B6.withchr_", winsize, ".cells_from_bin_analysis")
assertthat::assert_that(dir.exists(indir.raw))

# Load LDA objects --------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

infs.raw <- lapply(jmarks, function(jmark) list.files(indir.raw, pattern = paste0("B6-BM-", jmark, ".merged.NoCountThres.GeneTSS.Robj"), full.names = TRUE))

tm.result.lst <- lapply(infs, LoadGetTmResult)

inmain <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.dats <- lapply(jmarks, function(jmark) file.path(inmain, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))
dat.umap.long <- lapply(inf.dats, LoadUmap) %>%
  bind_rows() %>%
  dplyr::select(-repl, -techname)


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.umaps <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~mark) + 
  scale_color_manual(values = cbPalette)


# Load TSS ----------------------------------------------------------------

OldToNewName <- function(old.cname, i = 18){
  new.name <- paste(strsplit(old.cname, "\\.")[[1]][[i + 0]], strsplit(old.cname, "\\.")[[1]][[i + 1]], strsplit(old.cname, "\\.")[[1]][[i + 2]], strsplit(old.cname, "\\.")[[1]][[i + 3]], strsplit(old.cname, "\\.")[[1]][[i + 4]], sep = "-")
  return(new.name)
}
count.mat.tss <- lapply(infs.raw, function(x){
  load(x, v=T)
  # rename colnamescou
  colnames(count.dat$counts) <- sapply(names(colnames(count.dat$counts)), OldToNewName)
  # xout <- CollapseRowsByGene(count.dat$counts, as.long = TRUE)
  return(CollapseRowsByGene(count.dat$counts, as.long = FALSE))
})


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

dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat, stringsAsFactors = FALSE), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

# Find cell-type specific genes -------------------------------------------

# set up louvs
count.mats.louv <- lapply(count.mat.tss, function(count.mat) CollapseMatByLouvains(count.mat, dat.umap.long %>% mutate(louvmark = paste(mark, louvain, sep = "_"))))




# set up mats
jmark <- "H3K9me3"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
X <- count.mats.louv[[jmark]]
# normalize count size
ctype.metadata <- data.frame(ctype = colnames(X), mark = c(rep(jmark, ncol(X))))
rownames(ctype.metadata) <- ctype.metadata$ctype
dds <- DESeqDataSetFromMatrix(countData = X, colData = ctype.metadata, design = ~ 1)
vsd <- vst(dds)

X.norm <- normalize.quantiles(assay(vsd))
rownames(X.norm) <- rownames(assay(vsd)); colnames(X.norm) <- colnames(assay(vsd))
# X.norm <- log2(NormalizeMatrix(assay(vsd)))



print(unique(dat.norm.long$celltype))
jtop <- 500
from.top <- TRUE
# jctype <- "nucleate_erythrocyte"
# jctype <- "lymphocyte_of_B_lineage"
jctype <- "granulocyte"
# jctype <- "monocyte"

if (from.top){
  jsub <- subset(dat.norm.long, celltype == jctype) %>%
    group_by(celltype) %>%
    arrange(desc(zscore))
} else {
  jsub <- subset(dat.norm.long, celltype == jctype) %>%
    group_by(celltype) %>%
    arrange(zscore)
}
assertthat::assert_that(nrow(jsub) > 0)
jgenes <- jsub$gene[1:jtop]

# do H3K4me1 easy
X.norm.filt <- X.norm[rownames(X.norm) %in% jgenes, ] 
X.long <- data.table::melt(data.frame(gene = rownames(X.norm.filt), X.norm.filt)) %>%
  group_by(gene) %>%
  mutate(value = scale(value, center = TRUE, scale = TRUE))

m.box <- ggplot(X.long, aes(x = variable, y = value)) + geom_boxplot() + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle(jctype, paste("Ngenes:", nrow(X.norm.filt), "from.top", from.top))
  
m.umap <- ggplot(dat.umap.long %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_wrap(~mark) + 
  scale_color_manual(values = cbPalette)

multiplot(m.box, m.umap, cols = 2)

# Take top from mat, then ask what genes are in the bulk  -----------------

print(m.umap)
keep.top <- 250
jlouv <- 4
jlouv <- c(4, 7, 9)
jlouv <- 5
jlouv <- 3
jlouv <- 2
jlouv <- 6
jlouv <- 4
jlouv <- c(2, 8)
jlouv <- 3

jlouv <- 4


jlouv <- 7

jlouv <- 3

if (length(jlouv) == 1){
  louvname <- paste(jmark, jlouv, sep = "_")
} else {
  louvname <- paste(sapply(jlouv, function(x) paste(jmark, x, sep = "_")), sep = "|")
}

X.diff <- melt(gene = rownames(X.norm), X.norm)
x1 <- subset(X.diff, grepl(louvname, Var2)) %>% 
  dplyr::select(-Var2) %>%
  dplyr::rename(foreground = value,
                gene = Var1) 
x2 <- subset(X.diff, Var2 != louvname) %>% 
  group_by(Var1) %>%
  summarise(value = mean(value)) %>%
  dplyr::rename(background = value, 
                gene = Var1)
X.diff <- left_join(x1, x2, by = "gene") %>%
  mutate(fc = foreground - background) %>%
  # arrange(desc(fc))
  arrange(fc)

plot(X.diff$background, X.diff$foreground)
abline(a = 0, b = 1, col = 'blue', cex = 5)

jgenes <- X.diff$gene[1:keep.top]
jsub <- subset(dat.norm.long, gene %in% jgenes)

m.boxplot <- ggplot(jsub, aes(x = celltype, y = zscore)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(paste("Keeptop:", keep.top), paste("Louv:", louvname))

multiplot(m.boxplot, m.umap, cols = 2)






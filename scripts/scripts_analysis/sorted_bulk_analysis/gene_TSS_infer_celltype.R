# Jake Yeung
# Date of Creation: 2019-04-24
# File: ~/projects/scchic/scripts/scripts_analysis/sorted_bulk_analysis/gene_TSS_infer_celltype.R
# infer celltype from TSS

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(umap)
library(tidytext)


library(nnet)
library(msgl)
library(doParallel); library(foreach)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")


jmark <- "H3K4me1"

jfac <- 10^6
jpseudo <- 0

# Get trajs mixed ---------------------------------------------------------

trajs.mixed.out <- GetTrajMixed()
trajs.mixed <- trajs.mixed.out$trajs.mixed
dat.umap.mixed <- trajs.mixed.out$dat.umap.mixed


# Load bulk ---------------------------------------------------------------

inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat <- fread(inf.bulkdat, sep = "\t")

colnames(dat) <- gsub(" ", "_", colnames(dat))

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))

# Load UMAP ---------------------------------------------------------------

inf.umap <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.umap))
load(inf.umap, v=T)

# Load exprs  -------------------------------------------------------------

tssdist <- 50000
jdate <- "2019-04-22"
Kvec <- "50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, "_GeneTSS.Dedup.", jdate, ".", tssdist, ".Robj")
# inf.lda <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-50_GeneTSS.Dedup.2019-04-21.20000.Robj"

assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

out.lda <- out.lda[[length(out.lda)]]
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


# Get gene expression vector for each cell? -------------------------------

# can we correlate cells?
# library(multinom)
# jfit <- readRDS("/Users/yeung/data/scchic/robjs/multinom_fit.0.999.rds")

exprs.long <- exprs.long %>%
  mutate(exprs.log = log2(exprs * jfac + jpseudo)) %>%
  group_by(cell) %>%
  mutate(zscore = scale(exprs.log, center = TRUE, scale = TRUE))

# add trajectory info

m1 <- PlotXYWithColor(exprs.long %>% filter(gene == "S100a8"), xvar = "umap1", yvar = "umap2", cname = "zscore")
print(m1)


# Take 500 from each? -----------------------------------------------------

celltypes <- c("Kit_and_Sca1-positive_hematopoietic_stem_cell", "T_cell", "granulocyte", "nucleate_erythrocyte", "natural_killer_cell", "lymphocyte_of_B_lineage", "fetal_liver_hematopoietic_progenitor_cell", "megakaryocyte")
# celltypes <- unique(dat.sub$CellType)

tstart <- Sys.time()
keepn <- 50
genes.keep <- subset(top.peaks, rnk < keepn)$gene

# build classifier
dat.sub <- subset(dat.long %>% dplyr::rename(gene = Gene_Name), gene %in% genes.keep & CellType %in% celltypes)
dat.sub$CellTypeY <- relevel(as.factor(dat.sub$CellType), ref = "fetal_liver_hematopoietic_progenitor_cell")
dat.sub$gname <- paste("gene", dat.sub$gene, sep = "")

dat.mat <- tidyr::spread(dat.sub %>% ungroup() %>% dplyr::select(gname, CellTypeY, zscore), key = gname, value = zscore)

X <- dat.mat %>% dplyr::select(-CellTypeY)
# get Intercept
X <- as.matrix(cbind(data.frame(Intercept = rep(1, nrow(X)), X)))

# X <- model.matrix(~ gene * zscore, data = dat.sub)
Y <- dat.mat$CellTypeY

jfit.tss <- msgl::cv(x = X, classes = Y, standardize = FALSE, alpha = 1, lambda=0.1, use_parallel = FALSE, fold = 2)


save(jfit.tss, file = paste0("~/data/scchic/robjs/multinom_fit.tss", keepn, ".rds"))
print(Sys.time() - tstart)


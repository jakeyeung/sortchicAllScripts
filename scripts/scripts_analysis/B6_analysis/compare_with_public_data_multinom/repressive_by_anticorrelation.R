
rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)

library(preprocessCore)
library(here())

library(hash)

library(scchicFuncs)

setwd(here())
source("scripts/Rfunctions/PlotFunctions.R")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- "H3K4me3"
# jmark <- "H3K4me1"
# jmark <- "H3K9me3"
jmark <- "H3K27me3"


# Constants ---------------------------------------------------------------

use.zscore <- TRUE

# Load data ---------------------------------------------------------------

# inf <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.RData"
inf <- "/Users/yeung/data/scchic/robjs/B6_objs/terms_filt_H3K4me3_bin_TRUE_k_50.genomewide_nofilt.stringent_filter.withDist.RData"
load(inf, v=T)

inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)
out.objs.stringent <- out.objs
tm.result.stringent <- posterior(out.objs$out.lda)

inf.lda.all <- "/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
load(inf.lda.all, v=T)
out.objs$H3K4me3 <- out.objs.stringent
tm.result.lst <- lapply(out.objs, function(x) posterior(x$out.lda))
tm.result.lst[["H3K4me3"]] <- tm.result.stringent



inf.objs <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.objs, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent$H3K4me3


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

dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

if (use.zscore){
  dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>% as.data.frame() 
} else {
  dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-zscore), key = "celltype", value = "exprs") %>% as.data.frame() 
}
rownames(dat.norm.zscore.mat) <- dat.norm.zscore.mat$gene
dat.norm.zscore.mat$gene <- NULL


# keep only bins that are close to the TSS
terms.keep.dat <- out2.df.closest %>%
  group_by(gene) %>%
  filter(dist.to.tss == min(dist.to.tss))
terms.closest <- terms.keep.dat$region_coord

terms.sum <- terms.filt %>%
  group_by(gene) %>%
  dplyr::filter(rnk == min(rnk))

# term2gene <- hash(terms.keep.dat$region_coord, terms.keep.dat$gene)
# gene2term <- hash(terms.keep.dat$gene, terms.keep.dat$region_coord)

term2gene <- hash(terms.sum$term, terms.sum$gene)
gene2term <- hash(terms.sum$gene, terms.sum$term)

genes.keep <- rownames(dat.mat)



# Find correlated celltypes -----------------------------------------------


dat.pca <- prcomp(t(dat.mat), center = TRUE, scale. = TRUE)
dat.proj <- t(dat.mat) %*% dat.pca$rotation %*% diag(dat.pca$sdev)

plot(dat.proj[, 1], dat.proj[, 2])
text(dat.proj[, 1], dat.proj[, 2], labels = rownames(dat.proj))





# Set up raw dataa --------------------------------------------------------


inf.raw <- paste0("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_", jmark, "_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData")
load(inf.raw, v=T)
# handle rownames
rownames(count.dat$counts) <- paste("chr", rownames(count.dat$counts), sep = "")

all.cells <- colnames(count.dat$counts)
names(all.cells) <- all.cells

# take most interesting topics
top.nterms <- 150
topics.keep <- topics.sum$topic
# topics.keep <- c("Topic_40", "Topic_39", "Topic_27",   # monocyte
#                  "Topic_19", "Topic_8", "Topic_34", 
#                  "Topic_11", "Topic_6", "Topic_49", 
#                  "Topic_14", "Topic_1", "Topic_42")

terms.all <- unique(subset(terms.sum %>% arrange(desc(weight)), rnk <= top.nterms)$term)
genes.all <- unlist(sapply(terms.all, function(x) term2gene[[x]]))

# remove ribosomal genes
print(paste("N genes before:", length(genes.all)))
genes.all <- genes.all[!grepl("^Rp", x = genes.all)]
print(paste("N genes after:", length(genes.all)))



# Get differential expression  --------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# granu versus all
louv.in <- c(2)
louvs <- unique(dat.umap.long.trajs[[jmark]]$louvain)

for (louv.in in louvs){
  
  m.louv <- ggplot(dat.umap.long.trajs[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + geom_point() +
    scale_color_manual(values = cbPalette) +
    theme_bw () + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  cells.in <- subset(dat.umap.long.trajs[[jmark]], louvain == louv.in)$cell
  cells.out <- subset(dat.umap.long.trajs[[jmark]], louvain != louv.in)$cell
  
  # merge cells and check  differences 
  xvec.in <- Matrix::rowSums(count.dat$counts[, cells.in])
  xvec.out <- Matrix::rowSums(count.dat$counts[, cells.out])
  
  xvec.in <- log2(xvec.in / sum(xvec.in) * 10^6 + 1)
  xvec.out <- log2(xvec.out / sum(xvec.out) * 10^6 + 1)
  # xvec.in <- xvec.in / sum(xvec.in)
  # xvec.out <- xvec.out / sum(xvec.out)
  
  x.logFC <- sort(xvec.out - xvec.in, decreasing = TRUE)
  # x.logFC <- sort(xvec.out/xvec.in, decreasing = TRUE)
  # label top 100
  jterms.lab <- names(x.logFC)[1:100]
  jterms.lab.bottom <- names(sort(x.logFC, decreasing = FALSE))[1:100]
  xvec.in.filt <- xvec.in[c(jterms.lab, jterms.lab.bottom)]
  xvec.out.filt <- xvec.out[c(jterms.lab, jterms.lab.bottom)]
  jgenes.lab <- sapply(c(jterms.lab, jterms.lab.bottom), function(x){
    if (!is.null(term2gene[[x]])){
      term2gene[[x]]
    } else {
      x
    }
  })
  
  plot(xvec.in, xvec.out, pch = 20)
  points(xvec.in.filt, xvec.out.filt, pch = 20, cex = 3, col = 'blue')
  text(xvec.in.filt, xvec.out.filt, labels = jgenes.lab)
  abline(a = 0, b = 1, col = 'green', lwd = 3)
  
  
  # Plot top genes onto this  -----------------------------------------------
  
  jctypes <- unique(dat.norm.long$celltype)
  jctype <- "granulocyte"
  
  pdf(paste0("/Users/yeung/data/scchic/pdfs/celltyping_analysis_debugging/K27me3_louv_", paste(louv.in, collapse = "-"), ".pdf"), useDingbats = FALSE)
  
  print(m.louv)
  
  for (jctype in jctypes){
    dat.top <- dat.norm.long %>%
      rowwise() %>%
      mutate(is.celltype = celltype %in% jctype) %>%
      group_by(is.celltype, gene) %>%
      summarise(exprs = mean(exprs)) %>%
      group_by(gene) %>%
      summarise(log2FC = exprs[2] - exprs[1]) %>%
      arrange(desc(log2FC))
    
    # plot onto the  scatter
    jgenes.keep.tmp <- as.character(dat.top$gene[1:100])
    
    jterms.keep <- unique(subset(terms.filt, gene %in% jgenes.keep.tmp)$term)
    jgenes.keep <- sapply(jterms.keep, function(x){
      if (!is.null(term2gene[[x]])){
        return(term2gene[[x]])
      } else {
        return(x)
      }
    })
    
    # plot
    xvec.in.genefilt <- xvec.in[jterms.keep]
    xvec.out.genefilt <- xvec.out[jterms.keep]
    
    jterms.bg <- names(xvec.in)[which(!names(xvec.in) %in% jterms.keep)]
    
    n.genes.up.fg <- length(which(x.logFC[unique(jterms.keep)] >= 0))
    n.genes.up.bg <- length(which(x.logFC[unique(jterms.bg)] >= 0))
    
    n.genes.down.fg <- length(which(x.logFC[unique(jterms.keep)] < 0))
    n.genes.down.bg <- length(which(x.logFC[unique(jterms.bg)] < 0))
    
    X <- matrix(c(n.genes.up.fg, n.genes.up.bg, n.genes.down.fg, n.genes.down.bg), nrow = 2, byrow = FALSE)
    # X[1,1] <- 1000
    fish.out <- fisher.test(X, alternative = "greater")
    print(fish.out)
    # fish.pval <- fish.out$p.value
    
    plot(xvec.in, xvec.out, pch = 20, main = paste(jctype, "\nPval:", signif(fish.out$p.value, digits = 2), ", OR:", signif(fish.out$estimate, digits = 2)))
    points(xvec.in.genefilt, xvec.out.genefilt, pch = 20, cex = 3, col = 'blue')
    # text(xvec.in.genefilt, xvec.out.genefilt, labels = jgenes.lab, col = 'red')
    abline(a = 0, b = 1, col = 'green', lwd = 3)
    
  }
  dev.off()
  
}







# Jake Yeung
# Date of Creation: 2020-04-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/compare_pseudobulk_with_celltype_specific_genesets_custom.R
# In standard diff analysis, sometimes the mean levels matter (sometimes it doesn't?)



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Constants ---------------------------------------------------------------



# define for DE
jmean.min <- 4
fc.min <- 1.5
zscore.min <- 0
# define for seurat
padj.max <- 0.0001

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.downsample.compare_de_genes.custom"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("seurat_genes_on_chic_pseudobulk.padjmax_", padj.max, ".fcmin_", fc.min, ".pdf"))
outpdf2 <- file.path(outdir, paste0("de_genes_on_chic_pseudobulk.meanmin_", jmean.min, ".fcmin_", fc.min, ".zscoremin_", zscore.min, ".pdf"))
# outpdf3 <- file.path(outdir, paste0("seurat_genes_on_chic_pseudobulk.padjmax_", padj.max, ".fcmin_", fc.min, ".merged.pdf"))

# Load DE pseudobulk ------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_pseudobulk_scrnaseq_downsampled.2020-04-27.EosinophilsKeep.AllGenes.RData"
load(inf.de, v=T)


# Load seurat output ------------------------------------------------------

inf.seurat <- "/home/jyeung/data/from_rstudioserver/zebrafish.negbinom.2020-04-26/diff_exprs_Chloe_seurat.full.ctypefilt.rds"
# inf.seurat <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/from_data/zebrafish.poisson.2020-04-26/diff_exprs_Chloe_seurat.full.rds"

de.seurat <- readRDS(inf.seurat) %>%
  mutate(avg_log2FC = avg_logFC / log(2))

print(unique(de.seurat$cluster))

plot(density(de.seurat$avg_logFC))

ggplot(de.seurat, aes(x = avg_logFC, fill = cluster)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~cluster, ncol = 1) + xlim(c(-0.5, 0.5))

ggplot(de.seurat, aes(x = p_val_adj)) + 
  geom_density() + 
  facet_wrap(~cluster) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Define DE genes  --------------------------------------------------------



jctypes <- unique(as.character(pbulk.ctypefilt.long$pbulk))
# add NA
jctypes <- c(jctypes, "zOtherGenes")
# jctypes <- unique(as.character(pbulk.long$pbulk))
names(jctypes) <- jctypes

print(jctypes)

jctype <- "lymphocytes"
jsub <- subset(pbulk.ctypefilt.long, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
  arrange(desc(abs(log2fc)))


# do for all ctypes
de.genes.lst <- lapply(jctypes, function(jctype){
  if (jctype == "zOtherGenes"){
    jsub.tmp <- subset(pbulk.ctypefilt.long, log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
      # jsub <- subset(pbulk.long, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
      arrange(desc(abs(log2fc)))
    genes.tmp <- unique(jsub$gene)
    jsub <- subset(pbulk.ctypefilt.long, !gene %in% genes.tmp)
  } else {
    jsub <- subset(pbulk.ctypefilt.long, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
      # jsub <- subset(pbulk.long, pbulk == jctype & log2p1counts >= jmean.min & log2fc >= fc.min & zscore >= zscore.min) %>%
      arrange(desc(abs(log2fc)))
    genes.tmp <- jsub$gene
  }
  genes.ens <- sapply(as.character(genes.tmp), function(g) strsplit(g, "_")[[1]][[1]])
  return(genes.ens)
})

lapply(de.genes.lst, function(x) length(x))

de.genes.seurat.lst <- lapply(jctypes, function(jctype){
  if (jctype == "zOtherGenes"){
    genes.vec.tmp <- as.character(subset(de.seurat, p_val_adj <= padj.max & avg_log2FC >= fc.min)$gene)
    genes.vec <- unique(subset(de.seurat, !gene %in% genes.vec.tmp)$gene)
  } else {
    genes.vec <- subset(de.seurat, cluster == jctype & p_val_adj <= padj.max & avg_log2FC >= fc.min)$gene
  }
  ens.vec <- sapply(as.character(genes.vec), function(g) strsplit(g, "-")[[1]][[1]])
})

lapply(de.genes.seurat.lst, function(x) length(x))


# Plot top genes in the LDA and pseudobulks ----------------------------------------------

# load pseudobulks scchicseq

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jmark <- jmarks[[1]]


pbulk.long.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.pbulk.chic <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/rdata_robjs/WKM_scChICseq_pseudobulk_downsampled.ens.", jmark, ".RData")
  load(inf.pbulk.chic, v=T)
  return(pbulk.long)
})

# load LDA 
m.boxplots.logexprs <- lapply(jmarks, function(jmark){
  m <- ggplot(pbulk.long.lst[[jmark]], aes(x = pbulk, y = log2cuts)) + geom_boxplot() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
})
multiplot(m.boxplots.logexprs[[1]], m.boxplots.logexprs[[2]], m.boxplots.logexprs[[3]], cols = 3)

m.boxplots.ncuts <- lapply(jmarks, function(jmark){
  m <- ggplot(pbulk.long.lst[[jmark]], aes(x = pbulk, y = ncuts)) + geom_boxplot() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
})
multiplot(m.boxplots.ncuts[[1]], m.boxplots.ncuts[[2]], m.boxplots.ncuts[[3]], cols = 3)


m.boxplots.log2fc <- lapply(jmarks, function(jmark){
  m <- ggplot(pbulk.long.lst[[jmark]], aes(x = pbulk, y = log2FC)) + geom_boxplot() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
})
multiplot(m.boxplots.log2fc[[1]], m.boxplots.log2fc[[2]], m.boxplots.log2fc[[3]], cols = 3)

m.boxplots.zscore <- lapply(jmarks, function(jmark){
  m <- ggplot(pbulk.long.lst[[jmark]], aes(x = pbulk, y = log2zscore)) + geom_boxplot() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
})
multiplot(m.boxplots.zscore[[1]], m.boxplots.zscore[[2]], m.boxplots.zscore[[3]], cols = 3)

# get a gene set?
jgenes <- de.genes.seurat.lst$erythrocytes
jgenes <- de.genes.seurat.lst$lymphocytes
jgenes <- de.genes.seurat.lst$zOtherGenes
jgenes <- de.genes.seurat.lst$HSPCs
m.boxplots.zscore.filt <- lapply(jmarks, function(jmark){
  m <- ggplot(pbulk.long.lst[[jmark]] %>% filter(ens %in% jgenes), aes(x = pbulk, y = log2zscore)) + geom_boxplot() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark)
})
multiplot(m.boxplots.zscore.filt[[1]], m.boxplots.zscore.filt[[2]], m.boxplots.zscore.filt[[3]], cols = 3)


# Create one plot for the three marks?  -----------------------------------

ctypes.keep <- c("eryth", "HSC", "lymph", "monocyte")

# ctypes: eryth, HSC, lymph, monocytes
pbulk.k4me1 <- subset(pbulk.long.lst$H3K4me1) %>%
  ungroup() %>%
  mutate(mark = "H3K4me1") %>%
  filter(pbulk %in% ctypes.keep)

pbulk.k4me3 <- subset(pbulk.long.lst$H3K4me3) %>%
  ungroup() %>%
  mutate(mark = "H3K4me3") %>%
  filter(pbulk %in% ctypes.keep)

pbulk.k27me3 <- subset(pbulk.long.lst$H3K27me3) %>%
  ungroup() %>%
  mutate(mark = "H3K27me3",
         pbulk = gsub("eryth2", "eryth", pbulk),
         pbulk = gsub("HSC2", "HSC", pbulk)) %>%
  filter(pbulk %in% ctypes.keep)

# combine it all?
pbulk.merge <- bind_rows(pbulk.k4me1, pbulk.k4me3, pbulk.k27me3) %>%
  mutate(pbulk = gsub("monocyte", "granulocyte", pbulk))



# do for all genesets
de.ctypes.all <- names(de.genes.lst)
de.ctypes.filt <- c("erythrocytes", "granulocytes", "HSPCs", "lymphocytes", "zOtherGenes")


pbulk.merge.genesets <- lapply(de.ctypes.filt, function(de.ctype){
  ens.keep <- de.genes.lst[[de.ctype]]
  jsub <- subset(pbulk.merge, ens %in% ens.keep)
  ngenes <- length(unique(jsub$ens))
  jsub$geneset <- paste0(de.ctype, ",N=", ngenes)
  return(jsub)
})  %>%
  bind_rows()

m.merged.DE <- ggplot(pbulk.merge.genesets, aes(y = log2zscore, x = pbulk, fill = mark)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste("DE from pseudobulk")) + 
  geom_vline(xintercept = 0) +  geom_hline(yintercept = 0, linetype = "dotted") + 
  ylab("Zscore") + facet_wrap(~geneset)
print(m.merged.DE)



pbulk.merge.genesets.seurat <- lapply(de.ctypes.filt, function(de.ctype){
  ens.keep <- de.genes.seurat.lst[[de.ctype]]
  jsub <- subset(pbulk.merge, ens %in% ens.keep)
  ngenes <- length(unique(jsub$ens))
  jsub$geneset <- paste0(de.ctype, ",N=", ngenes)
  return(jsub)
})  %>%
  bind_rows()

m.merged.seurat <- ggplot(pbulk.merge.genesets.seurat, aes(y = log2zscore, x = pbulk, fill = mark)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste("DE from Seurat")) + 
  geom_vline(xintercept = 0) +  geom_hline(yintercept = 0, linetype = "dotted") + 
  ylab("Zscore") + facet_wrap(~geneset)
print(m.merged.seurat)

# plot K4me3 vs K27me3 
jsub.seurat <- pbulk.merge.genesets.seurat %>%
  reshape2::dcast(formula = ens + pbulk + geneset ~ mark, value.var = "ncuts")
jsub.de <- pbulk.merge.genesets %>%
  reshape2::dcast(formula = ens + pbulk + geneset ~ mark, value.var = "ncuts")

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m.corr.seurat.genesets <- ggplot(jsub.seurat %>% filter(!grepl("zOtherGenes", geneset)), aes(x = H3K4me3, y = H3K27me3, color = geneset)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~pbulk) + 
  scale_color_manual(values = cbPalette)

m.corr.seurat <- ggplot(jsub.seurat, aes(x = H3K4me3, y = H3K27me3)) + 
  geom_point(color = "black", alpha = 0.1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~pbulk) + 
  scale_color_manual(values = cbPalette)

m.dens.k4me3.seurat <- ggplot(jsub.seurat, aes(x = H3K4me3)) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m.denslog.k4me3.seurat <- ggplot(jsub.seurat, aes(x = H3K4me3)) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m.dens.k27me3.seurat <- ggplot(jsub.seurat, aes(x = H3K27me3)) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m.denslog.k27me3.seurat <- ggplot(jsub.seurat, aes(x = log10(H3K27me3))) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


m.corr.de.genesets <- ggplot(jsub.de %>% filter(!grepl("zOtherGenes", geneset)), aes(x = H3K4me3, y = H3K27me3, color = geneset)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~pbulk) + 
  scale_color_manual(values = cbPalette)

m.corr.de <- ggplot(jsub.de, aes(x = H3K4me3, y = H3K27me3)) + 
  geom_point(color = "black", alpha = 0.1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~pbulk) + 
  scale_color_manual(values = cbPalette)

m.dens.k4me3.de <- ggplot(jsub.de, aes(x = H3K4me3)) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m.denslog.k4me3.de <- ggplot(jsub.de, aes(x = log10(H3K4me3))) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m.dens.k27me3.de <- ggplot(jsub.de, aes(x = H3K27me3)) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m.denslog.k27me3.de <- ggplot(jsub.seurat, aes(x = log10(H3K27me3))) + 
  geom_density() + facet_wrap(~pbulk) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# Plot seurat -------------------------------------------------------------

pdf(outpdf, useDingbats = FALSE)
for (jctype in jctypes){
  ens.keep <- de.genes.seurat.lst[[jctype]]
  if (length(ens.keep) < 2){
    print(paste("Too few celltypes for", jctype, "skipping"))
  }
  m.lst <- lapply(jmarks, function(jmark){
    pbulk.long <- pbulk.long.lst[[jmark]]
    jsub <- subset(pbulk.long, ens %in% ens.keep)
    m <- ggplot(jsub, aes(x = log2zscore, group = pbulk, fill = pbulk)) + geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
      facet_wrap(~pbulk, ncol = 1) + ggtitle(paste("Seurat", jmark, jctype), paste("Ngenes:", length(ens.keep))) +
      geom_vline(xintercept = 0) + xlab("Zscore")
    return(m)
  })
  JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)
}
print(m.merged.seurat)
print(m.corr.seurat.genesets)
print(m.corr.seurat)
print(m.dens.k4me3.seurat)
print(m.dens.k27me3.seurat)
print(m.denslog.k4me3.seurat)
print(m.denslog.k27me3.seurat)
dev.off()


# Plot DE -----------------------------------------------------------------




pdf(outpdf2, useDingbats = FALSE)
for (jctype in jctypes){
  ens.keep <- de.genes.lst[[jctype]]
  m.lst <- lapply(jmarks, function(jmark){
    pbulk.long <- pbulk.long.lst[[jmark]]
    jsub <- subset(pbulk.long, ens %in% ens.keep)
    m <- ggplot(jsub, aes(x = log2zscore, group = pbulk, fill = pbulk)) + geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
      facet_wrap(~pbulk, ncol = 1) + ggtitle(paste("DE:", jmark, jctype), paste("Ngenes:", length(ens.keep))) +
      geom_vline(xintercept = 0) + xlab("Zscore")
    return(m)
  })
  JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)
}
print(m.merged.DE)
print(m.corr.de.genesets)
print(m.corr.de.genesets)
print(m.dens.k4me3.de)
print(m.dens.k27me3.de)
print(m.denslog.k4me3.de)
print(m.denslog.k27me3.de)
dev.off()

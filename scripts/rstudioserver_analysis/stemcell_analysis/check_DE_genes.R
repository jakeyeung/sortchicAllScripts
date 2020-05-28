# Jake Yeung
# Date of Creation: 2020-05-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/check_DE_genes.R
# We need a good set of DE genes. Check Giladi or the sorted celltypes

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(Seurat)

library(forcats)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.pdf"

exprs.max <- 3
logfc.min <- 2
logfc.max <- 1
exprs.thres.hk <- 5
exprs.thres.ne <- exprs.max * 0.5

pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)

# Load pseudobulks --------------------------------------------------------

m2c <- MarkerToCelltype()


WithVpreb <- TRUE

if (!WithVpreb){
  markers.keep <- c("Car1", "core", "Siglech", "Prg2", "Ccl5", "Prss34", "Cd74", "Fcrla", "Ltf")
} else {
  markers.keep <- c("Car1", "core", "Siglech", "Prg2", "Ccl5", "Prss34", "Cd74", "Fcrla", "Ltf", "Vpreb1")  # add Vpreb1?
}
markers.keep.again <- c("Erythroblast", "HSCs", "Bcell", "Neutrophil")

inf.mean <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
load(inf.mean, v=T)

dat.sum.norm.quantnorm <- dat.sum.norm.quantnorm[, markers.keep]

colnames(dat.sum.norm.quantnorm) <- sapply(colnames(dat.sum.norm.quantnorm), function(x) m2c[[x]])

print(head(dat.sum.norm.quantnorm))

dat.sum.norm.quantnorm <- dat.sum.norm.quantnorm[, markers.keep.again]

print(head(dat.sum.norm.quantnorm))

pbulk.exprs.long <- data.frame(gene = rownames(dat.sum.norm.quantnorm), dat.sum.norm.quantnorm, stringsAsFactors = FALSE) %>%
  reshape2::melt(data = ., id.vars = "gene", variable.name = "pbulk", value.name = "exprs") %>%
  group_by(gene) %>%
  mutate(logfc = exprs - mean(exprs), 
         zscore = logfc / sd(exprs))

pbulk.exprs.long$gene <- sapply(pbulk.exprs.long$gene, function(x) ifelse(grepl(";", x), strsplit(x, ";")[[1]][[1]], x))



# Load seurat outputs -----------------------------------------------------

# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.meth_poisson.Downsampled.FourCtypes.rds"
inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.meth_poisson.Downsampled.FourCtypes.WithVpreb_TRUE.rds"
assertthat::assert_that(file.exists(inf))

dat.de <- readRDS(inf)

pvalmax <- 0.001
logfcmin <- 0.1

jclsts <- as.character(unique(dat.de$cluster))
names(jclsts) <- jclsts

maxgenes <- Inf
de.genes.seurat <- lapply(jclsts, function(jclst){
  jsub <- subset(dat.de, cluster == jclst & p_val_adj <= pvalmax & avg_logFC >= logfcmin) %>%
    arrange(desc(avg_logFC))
  # get max genes
  gvec <- jsub$gene
  gvec.filt <- gvec[1:min(length(gvec), maxgenes)]
  # if ';' split 
  gvec.filt <- sapply(gvec.filt, function(x) ifelse(grepl(";", x), strsplit(x, ";")[[1]][[1]], x))
  return(gvec.filt)
})

jgenes.all <- sapply(dat.de$gene, function(g){
  gstrip <- ifelse(grepl(";", g), strsplit(g, ";")[[1]][[1]], g)
}, USE.NAMES = FALSE)

jens.all <- Gene2Ensembl(jgenes.all, return.original = TRUE)

g2e <- hash::hash(jgenes.all, jens.all)

de.ens.seurat <- sapply(de.genes.seurat, function(gvec){
  evec <- sapply(gvec, function(g) AssignHash(g, jhash = g2e, null.fill = g))
})

head(subset(dat.de, cluster == "Neutrophil") %>% arrange(desc(avg_logFC)))
head(subset(dat.de, cluster == "Erythroblast") %>% arrange(desc(avg_logFC)))
head(subset(dat.de, cluster == "HSCs") %>% arrange(desc(avg_logFC)))
head(subset(dat.de, cluster == "Bcell") %>% arrange(desc(avg_logFC)))

lapply(de.genes.seurat, length)

pbulk.exprs.long$ens <- sapply(pbulk.exprs.long$gene, function(g) AssignHash(x = g, jhash = g2e, null.fill = g))

# Label pseudobulks  ------------------------------------------------------



# Mean should matter ------------------------------------------------------

# plot 
m.lst <- lapply(jclsts, function(jclst){
  m <- pbulk.exprs.long %>% 
    filter(gene %in% de.genes.seurat[[jclst]]) %>%
    ggplot(., aes(x = pbulk, y = exprs)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(de.genes.seurat[[jclst]])))
  return(m)
})
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)

m.lst <- lapply(jclsts, function(jclst){
  m <- pbulk.exprs.long %>% 
    filter(gene %in% de.genes.seurat[[jclst]]) %>%
    ggplot(., aes(x = pbulk, y = logfc)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(de.genes.seurat[[jclst]])))
  return(m)
})
multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)


# Load sorted RNAseq ------------------------------------------------------



inf.bulkdat <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/E-MTAB-3079-query-results.fpkms.tsv"
dat.orig <- fread(inf.bulkdat, sep = "\t")
dat <- dat.orig

colnames(dat) <- gsub(" ", "_", colnames(dat))
print(colnames(dat))

# keep only certain colnames and rename to match above
cnames.keep <- c("Gene_ID", "Gene_Name", 
                 "Kit_and_Sca1-positive_hematopoietic_stem_cell", "granulocyte", 
                 "lymphocyte_of_B_lineage", "nucleate_erythrocyte")
cnames.keep.rename <- c("Gene_ID", "Gene_Name",
                        "HSCs", "Neutrophil", "Bcell", "Erythroblast")
dat <- dat[, ..cnames.keep]
colnames(dat) <- cnames.keep.rename

# print(dat)
# # replace NAs with 0s
# dat[is.na(dat)] <- 0
# print(dat)

dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
  group_by(Gene_ID) %>%
  mutate(FPKM = replace_na(FPKM, 0)) %>%
  # filter(all(!is.na(FPKM))) %>%
  group_by(Gene_ID, CellType) %>%
  summarise(FPKM = sum(FPKM)) %>%
  group_by(Gene_ID) %>%
  mutate(logFPKM = log2(FPKM + 1),
         logFC = logFPKM - mean(logFPKM),
         zscore = logFC / sd(logFPKM))

dat.mat <- tidyr::spread(dat.long %>%
                           ungroup() %>%
                           # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>%
                           mutate(gene = Gene_ID) %>%
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

dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "CellType", value = "exprs", -gene) %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
         logFC = scale(exprs, center = TRUE, scale = FALSE))

# normalize across samples?
ggplot(dat.norm.long, aes(x = CellType, y = exprs)) + geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.norm.long, aes(x = CellType, y = logFC)) + geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.norm.long, aes(x = CellType, y = zscore)) + geom_boxplot() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Try to find eryth-specific genes?  --------------------------------------



m.lst.sorted <- lapply(jclsts, function(jclst){
  jsub <- dat.norm.long  %>% 
    filter(gene %in% de.ens.seurat[[jclst]])
  m <- ggplot(jsub, aes(x = forcats::fct_reorder(.f = CellType, .x = exprs, .fun = median), y = exprs)) + 
    geom_boxplot() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(de.genes.seurat[[jclst]])))
  return(m)
})
multiplot(m.lst.sorted[[1]], m.lst.sorted[[2]], m.lst.sorted[[3]], m.lst.sorted[[4]], cols = 2)


m.lst.sorted <- lapply(jclsts, function(jclst){
  jsub <- dat.norm.long  %>% 
    filter(gene %in% de.ens.seurat[[jclst]])
  m <- ggplot(jsub, aes(x = forcats::fct_reorder(.f = CellType, .x = exprs, .fun = median), y = logFC)) + 
    geom_boxplot() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(de.genes.seurat[[jclst]])))
  return(m)
})
multiplot(m.lst.sorted[[1]], m.lst.sorted[[2]], m.lst.sorted[[3]], m.lst.sorted[[4]], cols = 2)


# # Get eryth specific genes in another way  --------------------------------
# 
# # high zscore, but mean of eryth should be high, and others low
# 
# eryth.dat.norm.long <- dat.norm.long %>% 
#   rowwise() %>%
#   mutate(is.eryth = CellType == "Erythroblast") %>%
#   group_by(gene, is.eryth) %>%
#   summarise(exprs = mean(exprs)) %>%
#   group_by(gene) %>%
#   mutate(logFC = exprs[[2]] - exprs[[1]],
#          exprs.mean = mean(exprs)) %>%
#   arrange(desc(logFC))
#   
# exprs.max <- 3
# logfc.min <- 2
# 
# eryth.de <- subset(eryth.dat.norm.long, exprs.mean <= exprs.max & logFC > logfc.min)
# 
# print(length(unique(eryth.de$gene)))
# # plot outputs
# 
# de.ens.seurat$ErythroblastStringent <- unique(as.character(eryth.de$gene))
# 
# jclsts <- c(jclsts, "ErythroblastStringent" = "ErythroblastStringent")


# Do for all  -------------------------------------------------------------

de.ens.sorted.stringent <- lapply(jclsts, function(jclst){
  print(jclst)
  jsub <- dat.norm.long %>% 
    rowwise() %>%
    mutate(is.ctype = CellType == jclst) %>%
    group_by(gene, is.ctype) %>%
    summarise(exprs = mean(exprs)) %>%
    group_by(gene) %>%
    mutate(logFC = exprs[[2]] - exprs[[1]],
           exprs.other = exprs[[1]]) %>%
    ungroup() %>%
    arrange(desc(logFC)) %>%
    filter(exprs.other <= exprs.max & logFC >= logfc.min)
  return(unique(jsub$gene))
})

lapply(de.ens.sorted.stringent, length)


lapply(de.ens.sorted.stringent, length)


# define house keeping?
jsub.hk <- dat.norm.long %>% 
  rowwise() %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs),
            exprs.sd = sd(exprs)) %>%
  arrange(exprs.mean) %>%
  filter(exprs.mean >= exprs.thres.hk & exprs.sd <= logfc.max)
genes.hk <- unique(jsub.hk$gene)

jsub.ne <- dat.norm.long %>% 
  rowwise() %>%
  group_by(gene) %>%
  summarise(exprs.mean = mean(exprs),
            exprs.sd = sd(exprs)) %>%
  arrange(desc(exprs.mean)) %>%
  filter(exprs.mean <= exprs.thres.ne & exprs.sd <= logfc.max)
genes.ne <- unique(jsub.ne$gene)

print(dim(jsub.ne))


de.ens.sorted.stringent$HighExprs <- genes.hk
de.ens.sorted.stringent$LowExprs <- genes.ne

# define not expressed

# add to jclsts
jclsts <- c(jclsts, "HighExprs" = "HighExprs", "LowExprs" = "LowExprs")

# Replot ------------------------------------------------------------------

# # Giladi
# m.lst <- lapply(jclsts, function(jclst){
#   print(jclst)
#   gvec <- de.ens.seurat[[jclst]]
#   jsub <- pbulk.exprs.long %>% 
#     filter(ens %in% gvec)
#   m <- jsub %>%
#     ggplot(., aes(x = pbulk, y = exprs)) + 
#     geom_boxplot() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(unique(jsub$ens))))
#   return(m)
# })
# multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], m.lst[[5]], cols = 2)
# 
# # Sorted
# m.lst2 <- lapply(jclsts, function(jclst){
#   print(jclst)
#   gvec <- de.ens.seurat[[jclst]]
#   jsub <- dat.norm.long %>% 
#     filter(gene %in% gvec)
#   m <- jsub %>%
#     ggplot(., aes(x = CellType, y = exprs)) + 
#     geom_boxplot() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(unique(jsub$gene))))
#   return(m)
# })
# multiplot(m.lst2[[1]], m.lst2[[2]], m.lst2[[3]], m.lst2[[4]], m.lst2[[5]], cols = 2)

# # look at the genes
# subset(pbulk.exprs.long, ens %in% de.ens.seurat[[jclsts]]) %>%
#   arrange(desc(logfc))

# # look at Hbb
# subset(pbulk.exprs.long, grepl("Sox6", gene))
# subset(eryth.dat.norm.long, grepl("ENSMUSG00000051910", gene))
# subset(eryth.de, grepl("ENSMUSG00000051910", gene))


# Check stringent in other genesets ---------------------------------------


m.stringent <- lapply(jclsts, function(jclst){
  gvec <- de.ens.sorted.stringent[[jclst]]
  jsub <- dat.norm.long %>% 
    filter(gene %in% gvec)
  m <- jsub %>%
    ggplot(., aes(x = CellType, y = exprs)) + 
    geom_boxplot() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste("Geneset:", jclst), paste("Ngenes:", length(unique(jsub$gene))))
  return(m)
})
multiplot(m.stringent[[1]], m.stringent[[2]], m.stringent[[3]], m.stringent[[4]], m.stringent[[5]], m.stringent[[6]], cols = 3)


# Check if these stringent eryth .de is already in Seurat -----------------

dat.de$ens <- sapply(dat.de$gene, function(g) AssignHash(x = g, jhash = g2e, null.fill = g))


eryth.dat.de <- subset(dat.de, ens %in% de.ens.seurat$ErythroblastStringent)  %>%
  filter(cluster == "Erythroblast")

eryth.dat.de.all <- subset(dat.de, ens %in% c(de.ens.seurat$ErythroblastStringent, de.ens.seurat$Erythroblast)) %>%
  filter(cluster == "Erythroblast")


# Save outputs ------------------------------------------------------------

# rename and save
dat.sorted.norm.long <- dat.norm.long
dat.giladi.norm.long <- pbulk.exprs.long
dat.de.giladi.seurat <- dat.de
# eryth.dat.sorted.norm.long <- eryth.dat.norm.long
g2e.hash <- g2e


head(dat.norm.long)

# de.ens.seurat

dev.off()
save(dat.sorted.norm.long, dat.giladi.norm.long, dat.de.giladi.seurat, g2e.hash, de.ens.seurat, de.ens.sorted.stringent, file = outf)


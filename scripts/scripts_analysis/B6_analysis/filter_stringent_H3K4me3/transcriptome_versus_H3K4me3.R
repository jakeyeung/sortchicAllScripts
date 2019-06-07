# Jake Yeung
# Date of Creation: 2019-06-06
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/transcriptome_versus_H3K4me3.R
# Compare with transcriptome 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(JFuncs)
library(topicmodels)
source("scripts/Rfunctions/PlotFunctions.R")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- "H3K4me3"
jmethod <- "pearson"
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
  mutate(logFPKM = log2(FPKM + 1),
         zscore = scale(logFPKM, center = TRUE, scale = TRUE))


# Can we correlate genome-wide clusters? ----------------------------------

# create cluster-specific transcriptomes? 
# imput.mat <- t(posterior(out.objs.stringent$out.lda)$topics %*% posterior(out.objs.stringent$out.lda)$terms)
imput.mat <- t(posterior(out.objs[[jmark]]$out.lda)$topics %*% posterior(out.objs[[jmark]]$out.lda)$terms)

# keep only bins that are close to the TSS
terms.keep.dat <- out2.df.closest %>%
  group_by(gene) %>%
  filter(dist.to.tss == min(dist.to.tss))
terms.keep <- terms.keep.dat$region_coord

imput.mat <- imput.mat[terms.keep, ]

terms2gene <- hash(terms.keep.dat$region_coord, terms.keep.dat$gene)


# rename rows to gene names?
rownames(imput.mat) <- sapply(rownames(imput.mat), function(x) terms2gene[[x]])

genes.keep <- rownames(imput.mat)

# merge across cells
louvains <- sort(unique(dat.umap.long.trajs[[jmark]]$louvain))

jscale <- 10^6
jpseudo <- 0
imput.merged.lst <- lapply(louvains, function(l){
  cells.keep <- subset(dat.umap.long.trajs[[jmark]], louvain %in% l)$cell
  cells.keep.i <- colnames(imput.mat) %in% cells.keep
  imput.sub <- imput.mat[, cells.keep.i]
  cell.exprs <- rowSums(imput.sub)
  out.dat <- data.frame(Gene_Name = names(cell.exprs), exprs = cell.exprs, cluster = l) %>%
    mutate(logexprs = log10(exprs * jscale + jpseudo))
  return(out.dat)
}) 
  
# Correlate expression with bulk transcriptome  ---------------------------

dat.sub <- dat.long %>%
  filter(Gene_Name %in% genes.keep) %>%
  group_by(Gene_Name) %>%
  filter(mean(logFPKM) > 2.5)

# genes.keep.exprs <- dat.sub$Gene_Name

# do only highly expressed genes??
plot(density(dat.sub$logFPKM))

cor.dat <- lapply(imput.merged.lst, function(imput.merged){
  tx.merged <- left_join(dat.sub, imput.merged)
  tx.cor <- tx.merged %>%
    group_by(CellType, cluster) %>%
    summarise(jcor = cor(logFPKM, logexprs, method = jmethod))
})  %>%
  bind_rows()

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "D3D3D3") 
cor.dat$jcol <- cbPalette[as.numeric(cor.dat$cluster)]

# plot individually and order properly 

print(cor.dat$jcor)
jylim <- c(min(0, floor(min(cor.dat$jcor) / 0.01) * 0.01), ceiling(max(cor.dat$jcor) / 0.01) * 0.01)

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/celltypes_top_genes_stringent/transcriptome_correlations_stringent_", jmark, "_method_", jmethod, ".pdf"), useDingbats = FALSE)
PlotXYWithColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jsize = 5)
m1 <- ggplot(cor.dat, aes(x = CellType, y = jcor, fill = jcol)) + geom_bar(stat = "identity") + facet_wrap(~cluster) + theme_bw() +
  scale_fill_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ylab(jmethod) + xlab("") + ylim(jylim)
print(m1)
for (l in louvains){
  jsub <- cor.dat %>% 
    filter(cluster == l)
  jsub <- OrderDecreasing(jsub, jfactor = "CellType", jval = "jcor")
  m1 <- ggplot(jsub, aes(x = CellType, y = jcor, fill = jcol)) + geom_bar(stat = "identity") + facet_wrap(~cluster) + theme_bw() +
    scale_fill_identity() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ylab(jmethod) + xlab("") + ylim(jylim)
  print(m1)
}
dev.off()

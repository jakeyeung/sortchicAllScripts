# Jake Yeung
# Date of Creation: 2020-12-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3b-make_umaps_bone_marrow_K27me3_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)


hubprefix <- "/home/jyeung/hub_oudenaarden"



# Load gene sets ----------------------------------------------------------


inf.genes.k4me1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt"
dat.genes.k4me1 <- fread(inf.genes.k4me1)

inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
dat.genes <- fread(inf.genes)


# Process genes -----------------------------------------------------------



dat.genes.all <- bind_rows(dat.genes, subset(dat.genes.k4me1, select = c(-topic)))

dat.genes.sub <- subset(dat.genes, !(jset == "Basophils" & rnk > 100))
dat.genes.sub <- subset(dat.genes.sub, !(jset == "pDCs" & rnk > 100))

dat.genes.sub <- subset(dat.genes, rnk < 150)

# swap basophils
dat.genes.sub <- subset(dat.genes.sub, jset != "Basophils" & jset != "pDCs")

basos.k4me3 <- subset(dat.genes, jset %in% c("Basophils"))$gene
basos.k4me1 <- subset(dat.genes.k4me1, jset %in% c("Basophils"))$gene

pdcs.k4me3 <- subset(dat.genes, jset %in% c("pDCs"))$gene
pdcs.k4me1 <- subset(dat.genes.k4me1, jset %in% c("pDCs"))$gene

basos.manual <- c("Il4", "Il6", "Cpa3", "Il1r1")
basos.rname <- subset(dat.genes.all, symbol %in% basos.manual)$gene
basos.both <- c(intersect(basos.k4me1, basos.k4me3), basos.rname)
pdcs.both <- intersect(pdcs.k4me1, pdcs.k4me3)

dat.to.add <- subset(dat.genes.all, gene %in% c(basos.both, pdcs.both))

dat.genes.sub.join <- bind_rows(dat.genes.sub, dat.to.add)

dat.genes.sub.join$jset <- factor(dat.genes.sub.join$jset, levels = c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs"))

dat.genes.sub.join <- dat.genes.sub.join %>%
  arrange(jset, rnk)



# Load Giladi  ------------------------------------------------------------

inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)

m2c <- MarkerToCelltype()
dat.public$celltype.annot <- sapply(dat.public$celltype, function(x) m2c[[x]])

# jgenes <- unique(subset(dat.genes.sub.join, jset == "Basophils")$symbol)
jgenes <- unique(subset(dat.genes.sub.join, jset == "")$symbol)
dat.sub <- subset(dat.public, gene %in% jgenes)

ggplot(dat.sub, aes(x = forcats::fct_reorder(.f = celltype, .x = zscore, .fun = median, .desc = TRUE), y = zscore)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  xlab("")



# Do heatmap  -------------------------------------------------------------

library(data.table)
mat.public <- subset(dat.public, gene %in% unique(dat.genes.sub.join$symbol)) %>%
  data.table::dcast(., formula = gene ~ celltype, value.var = "zscore")
rownames(mat.public) <- mat.public$gene
mat.public$gene <- NULL

pdf("/home/jyeung/hub_oudenaarden/jyeung/tmp/test.pdf")
heatmap3::heatmap3(as.matrix(mat.public[dat.genes.sub.join$symbol, ]), Rowv = NA, Colv = NA, scale = "row",  revC = TRUE)
dev.off()

# order columns appropriatelyu 
ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
ctypes.public <- c("Car1", "Fcrla", "Ccl5", "Ltf", "Prg2", "Siglech", "Cd74", "core")
ctypes.hash <- hash::hash(ctypes.public, ctypes)
ctypes2color <- hash::hash(dat.genes$jset, dat.genes$colorcode)

dat.public.sub <- subset(dat.public, celltype %in% ctypes.public) %>%
  rowwise() %>%
  mutate(jset = AssignHash(x = celltype, jhash = ctypes.hash, null.fill = x),
         colorcode = AssignHash(x = jset, jhash = ctypes2color, null.fill = jset))

mat.public.sub <- data.table::dcast(dat.public.sub, formula = gene ~ jset, value.var = "zscore")
rownames(mat.public.sub) <- mat.public.sub$gene
mat.public.sub$gene <- NULL



famous.genes.lst <- 
  list("Eryths" = c("Tal1", "Aqp1", "Gata1", "Comt", "Sphk1", "Hbb-bs", "Hbb-bt", "Sox6"),
       "Bcells" = c("Pax5", "Ebf1", "Cd79a", "Cd79b", "Snn", "Blnk", "Cd72", "Blk", "Kzf3", "Cd19"), 
       "NKs" = c("Stat4", "Tcf7", "Tbx21", "Cd3d", "Gimap4", "Gzma"), 
       "Granulocytes" = c("Cxcr2", "Mmp9", "Cd177", "Ncf4", "S100a9", "S100a8", "Ncf1", "Cebpe"),
       "Basophils" = c("Il4", "Il1r1", "Arap3", "Cpa3"), 
       "pDCs" = c("Irf8", "Selpg", "Itgb7", "Ccr9", "Unc93b1", "Siglech"), 
       "DCs" = c("Ccl9", "Apoe", "Nlrp3", "Csf1r"), 
       "HSPCs" = c("Hoxa9", "Hoxa7", "Hoxa3", "Meis1", "Runx2", "Kit", "Hlf", "Hoxa10", "Erg", "Cd34", "Hoxa6", "Gata3", "Hoxb4"))

famous.genes <- unlist(famous.genes.lst)




# pdf("/home/jyeung/hub_oudenaarden/jyeung/tmp/test.pdf")
pdf("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/compare_with_giladi/heatmap_geneset_from_chic_compare_giladi.pdf", useDingbats = FALSE)

jdedup <- dat.genes.sub.join %>%
  group_by(jset) %>%
  filter(!duplicated(symbol))

  indx <- which(jdedup$symbol %in% rownames(mat.public.sub))
  rkeep <- jdedup$symbol[indx]
  colskeep <- jdedup$colorcode[indx]
  mat.adj.tmp <- as.matrix(mat.public.sub[rkeep, ctypes])
  mat.adj.tmp <- t(apply(mat.adj.tmp, 1, function(jrow) DescTools::Winsorize(jrow, probs = c(0.01, 0.99))))
  mat.adj.tmp <- apply(mat.adj.tmp, 2, function(jcol) DescTools::Winsorize(jcol, probs = c(0.01, 0.99)))
  jlabrows <- rownames(mat.adj.tmp)
  # make NA if not in famous genes
  jlabrows[!jlabrows %in% famous.genes] <- NA
  heatmap3::heatmap3(mat.adj.tmp, Rowv = NA, Colv = NA, scale = "row", labRow = jlabrows, cexRow = 0.1, margin = c(5, 8), 
                     RowSideColors = colskeep, ColSideColors = sapply(ctypes, AssignHash, jhash = ctypes2color), 
                     main = "Gene sets defined by TSS ChIC, exprs from pseudobulk Giladi",
                     revC = TRUE)
  
dev.off()


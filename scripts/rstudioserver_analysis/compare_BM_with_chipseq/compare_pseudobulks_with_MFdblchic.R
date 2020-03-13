# Jake Yeung
# Date of Creation: 2020-03-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/compare_BM_with_chipseq/compare_pseudobulks_with_chipseq.R
# Cmopare pseudobulks with chipseq


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(GGally)
library(corrplot)
library(DescTools)
library(forcats)



# Load file ---------------------------------------------------------------

# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
jmark <- "H3K27me3"

indir <- paste0("/home/jyeung/hpc/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bigwigs_by_cluster.compare_with_MFdblchic.K27m3.again")
infs <- list.files(path = indir, pattern = "*.txt", full.names = TRUE)
infs.names <- sapply(infs, function(x){
  xsplit <- ClipLast(basename(x), jsep = "_")
  # xsplit <- strsplit(basename(x), split = "_")[[1]][[2]]  # /Users/yeung/data/dblchic/from_cluster/bigwigs/merged_bams.compare_with_Lara-Astiaso/H3K4me1_B_comparison.txt" -> B
}, USE.NAMES = FALSE)
names(infs) <- infs.names


# for (inf in infs){
cmats <- lapply(infs, function(inf){
  print(inf)
  dat <- fread(inf)
  # clean colnames
  colnames(dat) <- gsub("#", "", colnames(dat))
  colnames(dat) <- gsub("'", "", colnames(dat))
  colnames(dat) <- gsub(".bw", "", colnames(dat))
  colnames(dat) <- gsub(".sorted", "", colnames(dat))
  colnames(dat) <- gsub("_200119", "", colnames(dat))
  colnames(dat) <- gsub("all_BM", "allBM", colnames(dat))
  
  # Plot correlations? ------------------------------------------------------
  
  mat <- subset(dat, select = -c(chr, start, end))
  # remoev rows with NAs
  mat <- mat[complete.cases(mat), ]
  mat.win <- apply(mat, 2, function(jcol) Winsorize(jcol, probs = c(0, 0.98)))
  rownames(mat.win) <- rownames(mat)
  # ggpairs(data = mat)
  cmat <- cor(mat.win)
  corrplot(cmat, type = "upper")
  return(cmat[1, ])
})

cmats.dat <- do.call(rbind, cmats)
rownames(cmats.dat) <- names(cmats)

cmats.long <- data.frame(ctype.ref = rownames(cmats.dat), cmats.dat, stringsAsFactors = FALSE) %>%
  melt(., variable.name = "ctype.compare", value.name = "pcorr", id.vars = "ctype.ref") %>%
  filter(ctype.compare != "H3K4me1_B") %>%
  rowwise() %>%
  ungroup() %>%
  mutate(ctype.compare = as.factor(ctype.compare))

# rename things
# plot for each celltype a barplot of correlations
pdf(paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/correlations/correlate_", jmark, "_with_MFdblchic.pdf"), useDingbats = FALSE)
jctypes <- unique(cmats.long$ctype.ref)
for (jctype in jctypes){
  jsub <- cmats.long %>% filter(ctype.ref == jctype)
  m <- ggplot(jsub, aes(x = forcats::fct_reorder(.f = ctype.compare, .x = pcorr, .desc = TRUE), y = pcorr, fill = ctype.compare)) + 
    geom_col() + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
    ggtitle(paste(jmark, "reference:", jctype))
  print(m)
}
dev.off()  

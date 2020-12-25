# Jake Yeung
# Date of Creation: 2020-11-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/10-de_genes_k9me3.R
# 


rm(list=ls())

jstart <- Sys.time() 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)



library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

# Load UMAP  --------------------------------------------------------------




# Load differential exprs  ------------------------------------------------

# jmark <- "H3K9me3"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  print(jmark)
  
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_hiddendomains.", jmark, ".2020-11-17.newannot2.witherrors.RData")
  load(inf, v=T)
  
  ggplot(dat.annots.filt.mark, aes(x = umap1, y = umap2, color = Cluster)) + 
    geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  # Annotate domains to nearest gene  ---------------------------------------
  
  bins <- sapply(names(jfits.lst), function(x) strsplit(x, ";")[[1]][[1]])
  bins <- paste("chr", bins, sep = "")
  
  hubprefix <- "/home/jyeung/hub_oudenaarden"
  jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
  
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  
  bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = bins, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  
  dat.bins <- data.frame(bname = bins, rname = names(jfits.lst), stringsAsFactors = FALSE)
  
  plot(density(log10(abs(bins.annot$out2.df$dist.to.tss))))
  
  
  
  # Get differentially expressed regions  -----------------------------------
  
  
  jfit <- jfits.lst[[1]]
  
  jrow <- names(jfits.lst)[[1]]
  
  jrow.names <- names(jfits.lst); names(jrow.names) <- jrow.names
  
  dat.summary <- lapply(jrow.names, function(jrow){
    ctype.effects <- grep("^Cluster.*Estimate$", names(jfits.lst[[jrow]]), value = TRUE)
    ctype.stderrors <- grep("^Cluster.*StdError$", names(jfits.lst[[jrow]]), value = TRUE)
    fits.sub.effects <- unlist(jfits.lst[[jrow]][ctype.effects])
    fits.sub.stderrors <- unlist(jfits.lst[[jrow]][ctype.stderrors])
    dat.summary.tmp <- data.frame(param = names(fits.sub.effects), value = fits.sub.effects, se = fits.sub.stderrors, rname = jrow, stringsAsFactors = FALSE)
    return(dat.summary.tmp)
  }) %>%
    bind_rows()
  
  dat.pval <- lapply(jrow.names, function(jrow){
    pval <- jfits.lst[[jrow]]$pval
    dev.diff <- jfits.lst[[jrow]]$dev.diff
    dat.summary.tmp <- data.frame(pval = pval, dev.diff, rname = jrow, stringsAsFactors = FALSE)
    return(dat.summary.tmp)
  }) %>%
    bind_rows()
  
  dat.table <- lapply(jrow.names, function(jrow){
    jout <- data.frame(rname = jrow, as.data.frame(jfits.lst[[jrow]]), stringsAsFactors = FALSE)
    return(jout)
  }) %>%
    bind_rows()
  
  # Assign gene to rname  ---------------------------------------------------
  
  dat.bins.annot <- left_join(bins.annot$out2.df, dat.bins, by = c("region_coord" = "bname")) 
  dat.summary.annot <- left_join(dat.summary, subset(dat.bins.annot, select = c(rname, region_coord, tssname, gene, dist.to.tss))) %>%
    left_join(., dat.pval)  # add pval
  
  dat.table.annot <- left_join(dat.table, subset(dat.bins.annot, select = c(rname, region_coord, tssname, gene, dist.to.tss)))
  # remove duplicates
  dat.table.annot.dedup <- dat.table.annot %>%
    group_by(rname) %>%
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))
  
  
  ggplot(dat.summary, aes(x = value)) + 
    geom_density() 
  
  outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.downstream"
  fname <- paste0("poisson_fits_pval_and_effects_annotated.", jmark, ".", Sys.Date(), ".txt")
  fwrite(dat.table.annot.dedup, file = file.path(outdir, fname), na = "NA", sep = "\t", quote = FALSE)
  
  
}

print(Sys.time() - jstart)



# Get stem cell-specific effects ------------------------------------------

dat.summary.annot.filt <- subset(dat.summary.annot, abs(value) < 10) %>%
  group_by(rname) %>%
  filter(max(pval) < 1e-5)

head(dat.summary.annot %>% arrange(pval))

dat.summary.annot.filt.sum <- dat.summary.annot.filt %>%
  group_by(rname) %>%
  summarise(value.mean = mean(value), value.sd = sd(value), mlog10pval = mean(-log10(pval))) %>%
  arrange(desc(mlog10pval))


# Plot raw data -----------------------------------------------------------





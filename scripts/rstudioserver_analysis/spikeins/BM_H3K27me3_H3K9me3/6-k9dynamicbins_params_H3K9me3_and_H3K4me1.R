# Jake Yeung
# Date of Creation: 2021-01-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/6-k9dynamicbins_params_H3K9me3_and_H3K4me1.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load fits ---------------------------------------------------------------

jmark.ref <- "H3K9me3"

outs.lst <- lapply(jmarks, function(jmark){
  inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, pvals.long = pvals.long))
})

pvals.long.lst <- lapply(outs.lst, function(jout) jout$pvals.long)
params.long.lst <- lapply(outs.lst, function(jout) jout$params.long)


# Get celltype specific genes  --------------------------------------------

pval.cutoff <- 1e-10

bins.filt <- subset(pvals.long.lst[[jmark.ref]], pval < pval.cutoff)$bin

ggplot(params.long.lst[[jmark.ref]] %>% filter(bin %in% bins.filt), aes(x = estimate, fill = param)) + 
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  coord_cartesian(xlim = c(-5, 5)) + 
  theme_bw() + 
  facet_wrap(~param) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Compare H3K9me3 with other marks  ---------------------------------------

cnames.keep <- c("bin", "ClusterBcells.Estimate", "ClusterEryths.Estimate", "ClusterGranulocytes.Estimate")
cnames.keep2 <- cnames.keep[cnames.keep != "bin"]

# bins.common <- lapply(params.long.lst, function(jlong) unique(subset(jlong, bin %in% bins.filt)$bin))

bins.common <- Reduce(f = intersect, x = lapply(jmarks, function(jmark) unique(subset(params.long.lst[[jmark]], bin %in% bins.filt)$bin)))

bins.common.k4me1.k9me3 <- Reduce(f = intersect, x = lapply(jmarks[c("H3K4me1", "H3K9me3")], function(jmark) unique(subset(params.long.lst[[jmark]], bin %in% bins.filt)$bin)))

params.wide.lst <- lapply(jmarks, function(jmark){
  params.wide <- data.table::dcast(subset(params.long.lst[[jmark]], bin %in% bins.common), formula = "bin ~ param", value.var = "estimate")
  # assertthat::assert_that(all(cnames.keep %in% colnames(params.wide)))
  # params.wide.selected <- params.wide
  cnames.tmp <- colnames(params.wide)[colnames(params.wide) != "bin"]
  cnames.new <- c("bin", paste(cnames.tmp, jmark, sep = "."))
  colnames(params.wide) <- cnames.new
  return(params.wide)
})

params.wide.joined <- Reduce(f = left_join, x = params.wide.lst)


# Add pval ----------------------------------------------------------------


pvals.wide.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  pvals.sub <- subset(pvals.long.lst[[jmark]], bin %in% bins.common, select = c(bin, pval))
  cnames.new <- c("bin", paste0("pval.", jmark))
  colnames(pvals.sub) <- cnames.new
  return(pvals.sub)
})

pvals.wide.joined <- Reduce(f = left_join, x = pvals.wide.lst)

params.pvals.joined <- left_join(params.wide.joined, pvals.wide.joined, by = "bin")


# Make bed ----------------------------------------------------------------



dat.bins.filt <- data.frame(Chr = sapply(bins.filt, JFuncs::GetChromo),
                            Start = sapply(bins.filt, JFuncs::GetStart, returnAsInt = TRUE),
                            End = sapply(bins.filt, JFuncs::GetEnd, returnAsInt = TRUE), 
                            Name = bins.filt,
                            stringsAsFactors = FALSE) %>%
  rowwise() %>%
    mutate(Midpt = mean(c(Start, End)))

# fourth column needs to be gene, fifth column needs to be distance
bins.mid <- paste(dat.bins.filt$Chr, paste(dat.bins.filt$Midpt - 1, dat.bins.filt$Midpt + 1, sep = "-"), sep = ":")
bins.orig <- dat.bins.filt$Name

bins.hash <- hash::hash(bins.mid, bins.orig)

# Get closest gene -------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
bins.uniq <- unique(bins.mid)
dat.annot <- AnnotateCoordsFromList(coords.vec = bins.uniq, inf.tss = "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.100000.bed", txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

dat.annot2 <- dat.annot$regions.annotated

dat.annot2$bin <- sapply(bins.orig, function(x) AssignHash(x = x, jhash = bins.hash, null.fill = x))

params.pvals.joined.annot <- left_join(params.pvals.joined, dat.annot2)


# Write output ------------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/tables_H3K9me3_dynamic_bins"
outf.params <- file.path(outdir, paste0("log2FC_params_H3K9me3_dynamic_bins.closest_gene.", Sys.Date(), ".txt"))
outf.bed <- file.path(outdir, paste0("coords_H3K9me3_dynamic_bins.noname.", Sys.Date(), ".bed"))

fwrite(params.pvals.joined.annot, file = outf.params, sep = "\t")


# Write k9me3 bins  -------------------------------------------------------


fwrite(subset(dat.bins.filt, select = c(Chr, Start, End)), file = outf.bed, sep = "\t", col.names = FALSE)


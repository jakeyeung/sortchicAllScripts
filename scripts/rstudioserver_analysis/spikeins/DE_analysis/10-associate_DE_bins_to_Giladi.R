# Jake Yeung
# Date of Creation: 2021-03-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/DE_analysis/10-associate_DE_bins_to_Giladi.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(topGO)

library(CyclicGO)

# 
# 
# GetGOEnrichment <- function(genes.bg, genes.fg, fdr.cutoff, show.top.n = 8, ontology="BP", filter.GO.terms=FALSE){
#   # source(file.path(wd, "scripts/functions/AnalyzeGeneEnrichment.R"))
#   library(DBI)  # dbGetQuery() not found when loading topGO
#   library(topGO)
#   library(org.Mm.eg.db)
#   enrichment <- AnalyzeGeneEnrichment(genes.bg, genes.fg, FDR.cutoff = fdr.cutoff, which.ontology = ontology, return.GOdata = TRUE, filter.GO.terms = filter.GO.terms)
#   enrichment$minuslogpval <- -log10(as.numeric(enrichment$classicFisher))
#   
#   enrichment <- tryCatch({
#     enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
#   }, error = function(e) {
#     enrichment <- enrichment
#   })
#   
#   # enrichment <- OrderDecreasing(enrichment, jfactor = "Term", jval = "minuslogpval")
#   show.top.n.min <- min(nrow(enrichment), show.top.n)
#   if (show.top.n.min == 0) return(NULL)
#   enrichment <- enrichment[1:show.top.n.min, ]   # prevent taking more than you have enrichment
#   # unload packages
#   # sometimes topGO causes problems (unable to unload later), unload once you're done.
#   detach(name = "package:topGO", unload = TRUE)
#   detach(name = "package:org.Mm.eg.db", unload = TRUE)
#   return(enrichment)
# }



jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Filter DE ---------------------------------------------------------------



# filter DE bins 

infs.de.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

infs.high.lst <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned/High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dat.de.bins.lst <- lapply(infs.de.lst, function(jinf){
  fread(jinf)
})

dat.high.bins.lst <- lapply(infs.high.lst, function(jinf){
  fread(jinf)
})

de.bins.lst <- lapply(dat.de.bins.lst, function(jdat){
  jdat$CoordOriginal
})

high.bins.lst <- lapply(dat.high.bins.lst, function(jdat){
  jdat$CoordOriginal
})


# Load Giladi  -------------------------------------------------------------

inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
dat.public <- readRDS(inf.public)
dat.public$gene <- as.character(dat.public$gene)


# Separate HSPC high or low  ----------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jfits.lst.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})

params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>% 
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  # make params more readable
  params.dat.all$ctype <- params.dat.all$param
  params.dat.all$ctype <- gsub("Cluster", "", params.dat.all$ctype)
  params.dat.all$ctype <- gsub(".Estimate", "", params.dat.all$ctype)
  return(params.dat.all)
})



# Associate bins to gene --------------------------------------------------

# H3K27me3 



# jmarktmp <- "H3K27me3"
jmarktmp <- "H3K9me3"
# jmarktmp <- "H3K4me1"
# jmarktmp <- "H3K4me3"


jbins <- dat.de.bins.lst[[jmarktmp]]$CoordOriginal
jparams <- params.lst[[jmarktmp]]

jsum <- subset(jparams, bin %in% jbins) %>%
  group_by(bin) %>%
  summarise(estimate = mean(estimate))

plot(density(jsum$estimate))

jbins.lost <- subset(jsum, estimate < 0)$bin
jbins.gained <- subset(jsum, estimate > 0)$bin

maxdist <- 50000
jsub.lost <- subset(dat.de.bins.lst[[jmarktmp]], CoordOriginal %in% jbins.lost & abs(distanceToTSS) < maxdist)
jsub.gained <- subset(dat.de.bins.lst[[jmarktmp]], CoordOriginal %in% jbins.gained & abs(distanceToTSS) < maxdist)

jgene <- "Sulf1"
jgene <- "Msc"

jgenes.lost <- jsub.lost$SYMBOL
jgenes.random <- sample(unique(dat.public$gene), size = length(jgenes.lost))


jgenes.universe <- c(unique(dat.de.bins.lst[[jmarktmp]]$SYMBOL), unique(dat.high.bins.lst[[jmarktmp]]$SYMBOL))

# jgrep <- paste(jgenes, collapse = "|")

set.seed(0)
dat.public.sub <- subset(dat.public, gene %in% jgenes.lost) %>%
  mutate(type = paste0(jmarktmp, "_lost")) 

dat.public.random <- subset(dat.public, gene %in% jgenes.random)%>%
  mutate(type = paste0("random"))

ggplot(dat.public.sub, aes(x = celltype, y = exprs))  + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.public.random, aes(x = celltype, y = exprs))  + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.public.merged <- rbind(dat.public.sub, dat.public.random)

ggplot(dat.public.merged, aes(x = celltype, y = exprs, fill = type))  + 
  # geom_point() + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check gained  -----------------------------------------------------------


jgenes.gained <- jsub.gained$SYMBOL
jgenes.random.gained <- sample(unique(dat.public$gene), size = length(jgenes.gained))

# jgrep <- paste(jgenes, collapse = "|")

set.seed(0)
dat.public.sub.gained <- subset(dat.public, gene %in% jgenes.gained) %>%
  mutate(type = paste0(jmarktmp, "_gained")) 

dat.public.random.gained <- subset(dat.public, gene %in% jgenes.random.gained)%>%
  mutate(type = paste0("random"))

ggplot(dat.public.sub.gained, aes(x = celltype, y = exprs))  + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.public.random.gained, aes(x = celltype, y = exprs))  + 
  geom_point() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.public.merged.gained <- rbind(dat.public.sub.gained, dat.public.random.gained)

ggplot(dat.public.merged.gained, aes(x = celltype, y = exprs, fill = type))  + 
  # geom_point() + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Do GO term on K9me3 lost  -----------------------------------------------

go.lost <- GetGOEnrichment(genes.bg = jgenes.universe, genes.fg = jgenes.lost, fdr.cutoff = 0.05)
go.gained <- GetGOEnrichment(genes.bg = jgenes.universe, genes.fg = jgenes.gained, fdr.cutoff = 0.05)

# load Peter genes

inf.genes.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_peter/H3K27me3inHSClost.txt")
inf.genes.k9me3 <- file.path(hubprefix, "jyeung/data/scChiC/from_peter/H3K9me3inHSClost_log2of1.txt")

dat.genes.k27me3 <- fread(inf.genes.k27me3, header = FALSE)
dat.genes.k9me3 <- fread(inf.genes.k9me3, header = FALSE)

jgenes.universe.k27me3 <- unique(c(unique(dat.de.bins.lst[["H3K27me3"]]$SYMBOL), unique(dat.high.bins.lst[["H3K27me3"]]$SYMBOL)))
# jgenes.universe.k9me3 <- unique(c(unique(dat.de.bins.lst[["H3K9me3"]]$SYMBOL), unique(dat.high.bins.lst[["H3K9me3"]]$SYMBOL), unique(dat.genes.k9me3$V1)))
jgenes.universe.k9me3 <- c(unique(dat.de.bins.lst[["H3K9me3"]]$SYMBOL), unique(dat.high.bins.lst[["H3K9me3"]]$SYMBOL), unique(dat.genes.k9me3$V1))

go.check <- GetGOEnrichment(genes.bg = jgenes.universe.k9me3, genes.fg = dat.genes.k9me3$V1, )

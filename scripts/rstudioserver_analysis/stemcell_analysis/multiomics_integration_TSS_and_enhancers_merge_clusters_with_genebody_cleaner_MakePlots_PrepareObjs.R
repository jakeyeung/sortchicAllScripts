# Jake Yeung
# Date of Creation: 2020-04-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/multiomics_integration_TSS_and_enhancers_merge_clusters_with_genebody_cleaner.R
# description

# Clean up analysis after presentation
# Gene sets should be renamed to be meaningful
# Order the gene sets so things are on the diagonal 
# Compare density with a reference 
# absolute levels vs relative levels 

rm(list=ls())

jstart <- Sys.time()
library(reshape2)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DESeq2)

library(hash)

library(ggrepel)

m2c <- MarkerToCelltype()
m2c <- lapply(m2c, function(x) paste(x, "DEG", sep = "_"))
m2c <- lapply(m2c, function(x) gsub("NKcell", "NKcells", x))
m2c <- lapply(m2c, function(x) gsub("Bcell", "Bcells", x))

RenameLevels <- function(old.vec, conversion.lst, do.sort = TRUE){
  if (!do.sort){
    stop("Only sorting is implemented at the moment")
  }
  # order loadings in a sane way? 
  jlevels.old <- sort(as.character(unique(old.vec)), na.last = TRUE)
  jlevels.new <- sapply(jlevels.old, function(x) ifelse(!is.null(m2c[[x]]), m2c[[x]], NA))
  jlevels.hash <- hash::hash(jlevels.old[!is.na(jlevels.old)], jlevels.new[!is.na(jlevels.new)])
  # rename
  new.vec <- sapply(old.vec, function(x) jlevels.hash[[x]])
  new.vec <- factor(new.vec, levels = sort(jlevels.new))  # aphabeticica; 
  return(new.vec)
}

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#fa8c30", "#bd0a42", "#1347d5")
jdists <- c(1000L)
jsuffixs <- c(".celltypes_filt")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
# m2c <- MarkerToCelltype()

make.plots <- TRUE
overwrite <- FALSE

jdist <- 1000L
jsuffix <- jsuffixs[[1]]

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/proms_enhs_genebody_analysis2"
dir.create(outdir)

des.keep <- c("Car1", "Ccl5", "core", "Fcrla", "Ltf", NA)

inf.rdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/count_pseudobulks_enhs_proms_genebody.2020-04-02.RData"
outf.rdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/count_pseudobulks_enhs_proms_genebody.2020-04-02.ZscoreIsCtype.RData"
outf.rdata2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/count_pseudobulks_enhs_proms_genebody.2020-04-02.ZscoreIsCtype2.RData"

# Load objects ------------------------------------------------------------

load(inf.rdata, v=T)



# repare objets -----------------------------------------------------------


print(Sys.time() - jstart)
counts.pbulk.long.zscore.lst <- lapply(counts.pbulk.long.lst, function(jdat){
  jsub <- jdat %>%
    group_by(region_coord_full, de.ctype.choose) %>%
    mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
           logFC = scale(exprs, center = TRUE, scale = FALSE))
  # order loadings in a sane way? 
  jsub$de.ctype.choose <- RenameLevels(jsub$de.ctype.choose, m2c, do.sort = TRUE)
  jsub <- jsub %>%
    rowwise() %>%
    mutate(is.ctype = startsWith(x = as.character(de.ctype.choose), prefix = as.character(ctype))) %>%
    mutate(is.ctype = ifelse(is.na(is.ctype), FALSE, is.ctype))
  return(jsub)
})

print(Sys.time() - jstart)

save(counts.pbulk.long.lst, file = outf.rdata)

print(Sys.time() - jstart)
counts.pbulk.lst.zscore.allpseudos <- lapply(counts.pbulk.lst.allpseudos, function(jdat){
  jsub <- jdat %>%
    group_by(region_coord_full, de.ctype.choose) %>%
    mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
           logFC = scale(exprs, center = TRUE, scale = FALSE))
  # order loadings in a sane way? 
  jsub$de.ctype.choose <- RenameLevels(jsub$de.ctype.choose, m2c, do.sort = TRUE)
  jsub <- jsub %>%
    rowwise() %>%
    mutate(is.ctype = startsWith(x = as.character(de.ctype.choose), prefix = as.character(ctype))) %>%
    mutate(is.ctype = ifelse(is.na(is.ctype), FALSE, is.ctype))
  return(jsub)
})

print(Sys.time() - jstart)
save(counts.pbulk.long.lst, counts.pbulk.lst.allpseudos, file = outf.rdata2)

# 
# # Make plots --------------------------------------------------------------
# 
# # make density plots 
# 
# 
# 
# 
# # Check celltype specificity  ---------------------------------------------
# 
# 
# mshapes3 <- lapply(jmarks, function(jmark){
#   jsub <- subset(counts.pbulk.long.lst[[jmark]], de.ctype.choose %in% des.keep) %>%
#     group_by(region_coord_full) %>%
#     mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
#            logFC = scale(exprs, center = TRUE, scale = FALSE))
#   # order loadings in a sane way? 
#   jsub$de.ctype.choose <- RenameLevels(jsub$de.ctype.choose, m2c, do.sort = TRUE)
#   
#   m.exprs <- ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose, switch = "y") + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
#   m.zscore <- ggplot(jsub, aes(x = zscore, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose, switch = "y") + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist)) + 
#     geom_vline(xintercept = 0, linetype = "dotted", size = 0.5)
#   m.logFC <- ggplot(jsub, aes(x = logFC, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose, switch = "y") + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist)) + 
#     geom_vline(xintercept = 0, linetype = "dotted", size = 0.5)
#   return(list(m.exprs, m.zscore, m.logFC))
# })
# print(mshapes3)
# 
# 
# # Do some statistics comparing NA with others? ----------------------------
# 
# jlst <- list(c("H3K4me1", "enhancer"), c("H3K4me3", "promoter"), c("H3K27me3", "genebody"), c("H3K4me3", "enhancer"), c("H3K4me3", "genebody"))
# # jlst <- list(c("H3K4me3", ""))
# names(jlst) <- sapply(jlst, function(x) paste(x, collapse = "_"))
# 
# jsub.lst <- lapply(jlst, function(jvec){
#   jmark <- jvec[[1]]
#   jbtype <- jvec[[2]]
#   print(paste(jmark, jbtype))
#   
#   jsub <- subset(counts.pbulk.long.lst[[jmark]], de.ctype.choose %in% des.keep & biotype == jbtype) %>%
#     group_by(region_coord_full, mark, de.ctype.choose, biotype) %>%
#     mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
#            logFC = scale(exprs, center = TRUE, scale = FALSE)) %>%
#     ungroup() %>% 
#     mutate(de.ctype.choose = RenameLevels(de.ctype.choose, m2c))
#   
#   jsub <- jsub %>%
#     rowwise() %>%
#     mutate(is.ctype = startsWith(x = as.character(de.ctype.choose), prefix = as.character(ctype))) %>%
#     mutate(is.ctype = ifelse(is.na(is.ctype), FALSE, is.ctype))
#   # print(unique(jsub$is.ctype))
#   # jsub$is.ctype[is.na(jsub$is.ctype)] <- FALSE
#   return(jsub)
# })
# 
# m.compare <- lapply(names(jsub.lst), function(x){
#   jsub <- jsub.lst[[x]]
#   m <- ggplot(jsub, aes(x = logFC, fill = is.ctype, group = is.ctype)) + 
#     facet_wrap(~ctype) + geom_vline(xintercept = 0) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     geom_density(alpha = 0.5) + ggtitle(x)
# })
# print(m.compare)
# 
# print(unique(jsub$is.ctype))
# 
# m.promoter <- ggplot(jsub, aes(x = logFC, fill = is.ctype, group = is.ctype)) + 
#   facet_wrap(~ctype) + geom_vline(xintercept = 0) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_density(alpha = 0.5) + ggtitle(paste(jmark, jbtype))
# print(m.promoter)
# 
# m.enhancer <- ggplot(jsub, aes(x = logFC, fill = is.ctype, group = is.ctype)) + 
#   facet_wrap(~ctype) + geom_vline(xintercept = 0) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   geom_density(alpha = 0.5) + ggtitle(paste(jmark, jbtype))
# 
# 
# jsub$is.ctype <- 
# # model this
# fit.null <- lm(formula = logFC ~ 1 + is.ctype)
# 
# 
# x <- subset(jsub, ctype == "NKcells")$logFC
# x <- subset(jsub, ctype == "HSCs")$logFC
# y <- subset(jsub, ctype == "Bcells")$logFC
# ks.test(x, y)
# plot(density(x))
# plot(density(y))
# print(mean(x))
# print(mean(y))
# # put this in a model? 
# 
# 
# 
# mshapes3.nohsc <- lapply(jmarks, function(jmark){
#   jsub <- subset(counts.pbulk.long.lst[[jmark]], ctype != "HSCs" & de.ctype.choose %in% des.keep) %>%
#     group_by(region_coord_full) %>%
#     mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
#            logFC = scale(exprs, center = TRUE, scale = FALSE))
#   jsub$de.ctype.choose <- RenameLevels(jsub$de.ctype.choose, m2c, do.sort = TRUE)
#   
#   m.exprs <- ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
#   m.zscore <- ggplot(jsub, aes(x = zscore, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist)) + 
#     geom_vline(xintercept = 0, linetype = "dotted", size = 0.5)
#   m.logFC <- ggplot(jsub, aes(x = logFC, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist)) + 
#     geom_vline(xintercept = 0, linetype = "dotted", size = 0.5)
#   return(list(m.exprs, m.zscore, m.logFC))
# })
# print(mshapes3.nohsc)
# 
# 
# 
# 
# mshapes3.all <- lapply(jmarks, function(jmark){
#   jsub <- subset(counts.pbulk.lst.allpseudos[[jmark]], de.ctype.choose %in% des.keep) %>%
#     group_by(region_coord_full) %>%
#     mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
#            logFC = scale(exprs, center = TRUE, scale = FALSE))
#   jsub$de.ctype.choose <- RenameLevels(jsub$de.ctype.choose, m2c, do.sort = TRUE)
#   
#   m.exprs <- ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
#   m.zscore <- ggplot(jsub, aes(x = zscore, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist)) + 
#     geom_vline(xintercept = 0, linetype = "dotted")
#   m.logFC <- ggplot(jsub, aes(x = logFC, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
#     facet_grid(ctype ~ de.ctype.choose) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist)) + 
#     geom_vline(xintercept = 0, linetype = "dotted")
#   return(list(m.exprs, m.zscore, m.logFC))
# })
# print(mshapes3.all)



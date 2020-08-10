# Jake Yeung
# Date of Creation: 2020-07-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/check_lm_var_correction_output.with_2D_plot.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)
library(ggrepel)

FitClstVar <- function(jrow, jmeta){
  fit.input <- data.frame(exprs = jrow, xvar = jmeta$xvar, clst = jmeta$clst, stringsAsFactors = TRUE)
  jfit <- lm(exprs ~ xvar:clst + clst, fit.input)
  return(jfit)
}

# FitClstVar.null <- function(jrow, jmeta){
#   fit.input <- data.frame(exprs = jrow, xvar = jmeta$xvar, clst = jmeta$clst, stringsAsFactors = TRUE)
#   jfit <- lm(exprs ~ xvar:clst + 1, fit.input)
#   return(jfit)
# }


make.plots <- TRUE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/pdfs_from_WKM_BM_merged"

outf <- file.path(outdir, paste0("plot_varcorrection_and_logfc_examples.", Sys.Date(), ".2D_and_volcanos.fixed3.pdf"))

if (make.plots){
  pdf(file = outf, useDingbats = FALSE)
}

# Load gene sets ----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

# Load inputs -------------------------------------------------------------

# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits.2020-06-25.RData"
# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits.2020-06-26.CleanUpEryth.RData"
inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits.2020-07-14.WithNull.RedoPval.RData"
load(inf, v=T)

dat.params.all.lst <- lapply(dat.params.all.lst, function(jdat){
  jsub.int <- subset(jdat, params.mean == "clstHSC")
  jsub.int.tomerge <- data.frame(rname = jsub.int$rname, intercept = jsub.int$jmean)
  jdat.merged <- left_join(jdat, jsub.int.tomerge)
})

# Check DE of gene sets  --------------------------------------------------

dat.params.all.annot <- lapply(jmarks, function(jmark){
  jdat <- dat.params.all.lst[[jmark]]
  jdat$gene <- sapply(jdat$rname, function(x) strsplit(x, ";")[[1]][[2]])
  jdat$ens <- sapply(jdat$gene, AssignHash, g2e.hash)
  jdat$mark <- jmark
  jdat <- subset(jdat, !is.na(ens))
  return(jdat)
}) %>%
  bind_rows()


# do by gene lists 
jens.vec.names <- names(de.ens.sorted.stringent)
names(jens.vec.names) <- jens.vec.names

params.bygsets <- lapply(jens.vec.names, function(jname){
  ens.vec <- de.ens.sorted.stringent[[jname]]
  jsub <- subset(dat.params.all.annot, ens %in% ens.vec)
  jsub$gset <- jname
  jsub <- jsub %>%
    filter(params.mean != "clstHSC") %>%
    mutate(logfc = jmean - intercept)
  return(jsub)
}) %>%
  bind_rows()



cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(params.bygsets, aes(x = logfc, fill = gset)) + facet_grid(mark~params.mean) + geom_density(alpha = 0.25) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dotted") + scale_fill_manual(values = cbPalette)

for (jmark in jmarks){
  m <- ggplot(params.bygsets %>% filter(mark == jmark), aes(x = logfc, fill = gset)) + facet_grid(gset~params.mean) + geom_density(alpha = 0.25) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_vline(xintercept = 0, linetype = "dotted") + scale_fill_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}

# plot celltype specific
jclsts <- sort(unique(params.bygsets$params.mean)); names(jclsts) <- jclsts
print(unique(params.bygsets$gset))
jgsets <- c("Erythroblast", "Neutrophil", "Bcell"); names(jgsets) <- jclsts

for (jmark in jmarks){
  print(jmark)
  for (jclst in jclsts){
    print(jclst)
    other.clst <- paste(jclsts[which(jclsts != jclst)], collapse = ",")
    jgset <- jgsets[[jclst]]
    params.bygsets.tmp <- params.bygsets %>%
      filter(gset == jgset & mark == jmark) %>%
      rowwise() %>%
      mutate(params.mean = ifelse(params.mean == jclst, jclst, other.clst))
    params.bygsets.tmp$params.mean <- factor(params.bygsets.tmp$params.mean, levels = c(jclst, other.clst))
    m.clst <- ggplot(params.bygsets.tmp, mapping = aes(x = logfc, fill = params.mean)) + 
      geom_density(alpha = 0.33) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      ggtitle(paste(jmark, jclst, jgset)) 
    print(m.clst)
  }
}
  



# Plot some examples ------------------------------------------------------


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")




# Neutrophil: S100a7a? 
jgene <- "S100a7a"

jgene <- "Retnlg"


jgene <- "S100a7a"

jmark <- "H3K4me1"
jmark <- "H3K27me3"

jmark <- "H3K4me3"

jgene <- "Irf4"
jgene <- "Tal1"
jgene <- "S100a7a"



jgene <- "Cebpd"

jgene <- "Pax5"

jgene <- "Irf4"

jgene <- "Ebf1"

jgene <- "Sox6"

jgene <- "Hbb-y"


jgene <- "Hoxa4"
jgene <- "S100a8"

m.all <- ggplot(dat.params.all.annot %>% filter(params.mean != "clstHSC"), mapping = aes(x = jmean - intercept, fill = mark)) + geom_density()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + geom_vline(xintercept = 0) + xlab("logFC")
print(m.all)

m.all.byclst <- ggplot(dat.params.all.annot %>% filter(params.mean != "clstHSC"), mapping = aes(x = jmean - intercept, fill = params.mean)) + geom_density(alpha = 0.33)  +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + geom_vline(xintercept = 0) + xlab("logFC")
print(m.all.byclst)

m.all.byclst.facet <- ggplot(dat.params.all.annot %>% filter(params.mean != "clstHSC"), mapping = aes(x = jmean - intercept, fill = params.mean)) + geom_density(alpha = 0.33)  +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(params.mean~mark) + geom_vline(xintercept = 0) + xlab("logFC")
print(m.all.byclst.facet)




# Make 2D plots -----------------------------------------------------------

# outpdf <- "~/data/rsessions/check_lm_var_correction.with_2D_plots.genelab.WithVolcano.pdf"
# pdf(outpdf, useDingbats = FALSE)

params.bygsets.nodupes <- params.bygsets %>%
  group_by(ens, mark) %>%
  filter(length(unique(gset)) == 1)
  
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ctypes.end <- list("clstgranu", "clstlymph", "clsteryth"); names(ctypes.end) <- ctypes.end
gset.specs <- list("Neutrophil", "Bcell", "Erythroblast", "HighExprs"); names(gset.specs) <- ctypes.end
gset.others <- list(c("Bcell", "Erythroblast"), c("Erythroblast", "Neutrophil"), c("Bcell", "Neutrophil"), c("LowExprs")); names(gset.others) <- ctypes.end
gset.others.names <- lapply(gset.others, function(jgset) paste(jgset, collapse = "&")); names(gset.others.names) <- ctypes.end
# gset.HSCs <- "HSCs"
gset.HSCs <- list(c("HSCs"), c("HSCs"), c("HSCs"), c("HSCs")); names(gset.HSCs) <- ctypes.end


jmarks.keep <- c("H3K4me3", "H3K27me3")
print(unique(dat.params.all.annot$params.mean))


ref.clst <- "clstHSC"
# jctype <- "clstHSC"
jctype <- "clstgranu"
jctype <- "clsteryth"
jctype <- "clstlymph"

jctypes <- c("clstgranu", "clsteryth", "clstlymph"); names(jctypes) <- jctypes

for (jctype in jctypes){
  
  
  jtitle <- paste0(ref.clst, " -> ", jctype)
  print(jtitle)
  # jgset.name <- "Erythroblast"
  # jgset <- subset(params.bygsets, gset == jgset.name)$ens
  
  # plot Spec vs NotSpec genesets
  jsub.wide.logfc.Spec <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.specs[[jctype]]) %>%
    reshape2::dcast(data = ., formula = "gene + ens ~ mark", value.var = "logfc")  %>%
    mutate(gset = gset.specs[[jctype]])
  jsub.wide.logfc.NonSpec <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.others[[jctype]]) %>% 
    reshape2::dcast(data = ., formula = "gene + ens ~ mark", value.var = "logfc")  %>%
    mutate(gset = gset.others.names[[jctype]])
  
  jsub.wide.logfc.SpecNonSpec <- bind_rows(jsub.wide.logfc.Spec, jsub.wide.logfc.NonSpec)
  jsub.wide.logfc.SpecNonSpec$gset <- factor(jsub.wide.logfc.SpecNonSpec$gset, levels = c(gset.specs[[jctype]], gset.others.names[[jctype]]))
  
  # label gene if in top 10 of log2fc
  jsub.wide.logfc.SpecNonSpec <- jsub.wide.logfc.SpecNonSpec %>%
    rowwise() %>%
    mutate(veclength = sqrt(H3K4me3^2 + H3K27me3^2)) %>%
    group_by(gset) %>%
    mutate(gene.lab = ifelse(veclength >= quantile(veclength, probs = 0.95), gene, NA))
  
  
  m <- ggplot(jsub.wide.logfc.SpecNonSpec, aes(x = H3K4me3, y = H3K27me3, color = gset, label = gene.lab)) + 
    geom_point(alpha = 0.25)  + 
    geom_text_repel(color = 'black', size = 3, color = "grey60", alpha = 0.8) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    facet_wrap(~gset) + ggtitle(jtitle)
  print(m)
  
  
  # plot Spec vs HSC genesets
  jsub.wide.logfc.Spec <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.specs[[jctype]]) %>%
    reshape2::dcast(data = ., formula = "gene + ens ~ mark", value.var = "logfc")  %>%
    mutate(gset = gset.specs[[jctype]])
  jsub.wide.logfc.HSC <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.HSCs[[jctype]]) %>% 
    reshape2::dcast(data = ., formula = "gene + ens ~ mark", value.var = "logfc")  %>%
    mutate(gset = gset.HSCs[[jctype]])
  
  jsub.wide.logfc.SpecHSC <- bind_rows(jsub.wide.logfc.Spec, jsub.wide.logfc.HSC)
  jsub.wide.logfc.SpecHSC$gset <- factor(jsub.wide.logfc.SpecHSC$gset, levels = c(gset.specs[[jctype]], gset.HSCs[[jctype]]))
  
  # label gene if in top 10 of log2fc
  jsub.wide.logfc.SpecHSC <- jsub.wide.logfc.SpecHSC %>%
    rowwise() %>%
    mutate(veclength = sqrt(H3K4me3^2 + H3K27me3^2)) %>%
    group_by(gset) %>%
    mutate(gene.lab = ifelse(veclength >= quantile(veclength, probs = 0.95), gene, NA))
  
  
  m <- ggplot(jsub.wide.logfc.SpecHSC, aes(x = H3K4me3, y = H3K27me3, color = gset, label = gene.lab)) + 
    geom_point(alpha = 0.25)  +  
    geom_text_repel(color = 'black', size = 3, color = "grey60", alpha = 0.8) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    facet_wrap(~gset) + ggtitle(jtitle)
  print(m)
  
  # plot all genesets
  jgset.names <- unique(params.bygsets$gset)
  names(jgset.names) <- jgset.names
  
  jsub.wide.logfc.all <- lapply(jgset.names, function(jgset.name){
    dat <- subset(params.bygsets, mark %in% jmarks.keep & params.mean == jctype & gset == jgset.name) %>%
      reshape2::dcast(data = ., formula = "gene + ens ~ mark", value.var = "logfc") 
    dat$gset <- jgset.name
    return(dat)
  }) %>%
    bind_rows()
  
  m <- ggplot(jsub.wide.logfc.all, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    geom_point(alpha = 0.2)  + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_vline(xintercept = 0, alpha = 0.5, linetype = "dotted") + 
    geom_hline(yintercept = 0, alpha = 0.5, linetype = "dotted") + 
    scale_color_manual(values = cbPalette)  + 
    ggtitle(jtitle)
  print(m)
  
  
  # Add  2D arrow -----------------------------------------------------------
  
  ref.clst <- "clstHSC"
  # end.clst <- "clsteryth"
  end.clst <- jctype
  jgset1 <- gset.specs[[end.clst]]
  jclsts <- c(ref.clst, end.clst)
  names(jclsts) <- jclsts
  
  # jgset2 <- gset.others[[end.clst]]
  # jgset2 <- gset.HSCs[[end.clst]]
  
  jsub.wide.2d <- lapply(jclsts, function(jclst){
    jsub.wide <- subset(dat.params.all.annot, mark %in% jmarks.keep & params.mean == jclst) %>%
      reshape2::dcast(data = ., formula = "gene + ens ~ mark", value.var = "jmean")
    cnames.old <- grep(pattern = "H3K", colnames(jsub.wide), value = TRUE)
    cnames.new <- paste(cnames.old, ifelse(jclst == ref.clst, "RefClst", "EndClst"), sep = "_")
    colnames(jsub.wide) <- c("gene", "ens", cnames.new)
    return(jsub.wide)
  }) %>%
    Reduce(f = left_join, x = .) %>%
    rowwise() %>%
    mutate(H3K4me3_logFC = H3K4me3_EndClst - H3K4me3_RefClst, 
           H3K27me3_logFC = H3K27me3_EndClst - H3K27me3_RefClst)
  
  jfactor <- 0.3
  for (jgset2 in c(gset.others[[end.clst]], gset.HSCs[[end.clst]])){
    
    jgset.lst <- list(g1 = jgset1, g2 = jgset2)
    jgset.lst.names <- lapply(jgset.lst, function(x) paste(x, collapse = "&"))
    names(jgset.lst) <- jgset.lst.names
    jnames <- names(jgset.lst); names(jnames) <- jnames
    
    jsub.wide.2d.gset <- lapply(jnames, function(jname){
      jgset <- jgset.lst[[jname]]
      print(jgset)
      jsub <- jsub.wide.2d %>% filter(ens %in% subset(params.bygsets, gset %in% jgset)$ens)
      jsub$gset <- jname
      return(jsub)
    }) %>%
      bind_rows()
    jsub.wide.2d.gset$gset <- factor(jsub.wide.2d.gset$gset, levels = names(jgset.lst))
    
    m <- ggplot(jsub.wide.2d.gset, aes(x = H3K4me3_RefClst, xend = H3K4me3_RefClst + H3K4me3_logFC * jfactor, y = H3K27me3_RefClst, yend = H3K27me3_RefClst + H3K27me3_logFC * jfactor, color = gset)) + 
      geom_segment(arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.8, size = 0.1) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      xlab("H3K4me3 log signal") + 
      ylab("H3K27me3 log signal") + 
      facet_wrap(~gset) + ggtitle(jtitle)
    print(m)
  }
  
  
}






# Add volcanos ------------------------------------------------------------

head(params.bygsets)

# genome wide for each jtype
for (jctype in jctypes){
  jsub <- params.bygsets %>% filter(params.mean == jctype)
  jrange <- max(abs(jsub$logfc))
  jtitle <- paste0("HSPCs -> ", jctype)
  m <- ggplot(jsub, aes(x = logfc, y = -log10(pval))) + geom_point(alpha = 0.1)  + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_grid(mark ~ gset) + xlim(c(-jrange, jrange)) + 
    ggtitle(jtitle)
  print(m)
  m.dens <- ggplot(jsub, aes(x = logfc, fill = gset)) + geom_density()  + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_grid(mark ~ gset) + xlim(c(-jrange, jrange)) + 
    ggtitle(jtitle) + geom_vline(xintercept = 0, linetype = "dotted")
  print(m.dens)
  
  jsub.long.logfc.Spec <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.specs[[jctype]]) %>%
    mutate(gset = gset.specs[[jctype]])
  jsub.long.logfc.NonSpec <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.others[[jctype]]) %>% 
    mutate(gset = gset.others.names[[jctype]])
  jsub.long.logfc.HSCs <- subset(params.bygsets.nodupes, mark %in% jmarks.keep & params.mean == jctype & gset %in% gset.HSCs[[jctype]]) %>% 
    mutate(gset = gset.HSCs[[jctype]])
  jsub.long.logfc.merge <- bind_rows(jsub.long.logfc.Spec, jsub.long.logfc.NonSpec)
  jsub.long.logfc.merge.hscs <- bind_rows(jsub.long.logfc.Spec, jsub.long.logfc.HSCs)
  
  jsub.long.logfc.merge$gset <- factor(jsub.long.logfc.merge$gset, levels = c(gset.specs[[jctype]], gset.others.names[[jctype]]))
  jsub.long.logfc.merge.hscs$gset <- factor(jsub.long.logfc.merge.hscs$gset, levels = c(gset.specs[[jctype]], gset.HSCs[[jctype]]))
  
  m.vol <- ggplot(jsub.long.logfc.merge, aes(x = logfc, y = -log10(pval), color = gset)) + geom_point() + facet_grid(mark ~ gset) + 
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    ggtitle(jtitle) + xlim(c(-jrange, jrange))
  print(m.vol)
  m.vol.hscs <- ggplot(jsub.long.logfc.merge.hscs, aes(x = logfc, y = -log10(pval), color = gset)) + geom_point() + facet_grid(mark ~ gset) + 
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    ggtitle(jtitle) + xlim(c(-jrange, jrange))
  print(m.vol.hscs)
}

# dev.off()

# 
# ggplot(params.bygsets, aes(x = logfc)) + geom_density()  + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   facet_grid(mark ~ gset) + geom_vline(xintercept = 0)


# save.image(file = "~/data/rsessions/lm_var_correction_downstream_again.with_pval.RData")


# plot individual genes ---------------------------------------------------

# m.test <- ggplot(subset(params.bygsets.nodupes, mark == "H3K4me1"), aes(x = logfc, y = -log10(pval))) + facet_wrap(~gset) + geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# print(m.test)
# 
# jsub <- subset(params.bygsets.nodupes, mark == "H3K4me1" & abs(logfc) < 0.1 & -log10(pval) > 300) %>%
#   arrange(abs(logfc))
#   # arrange(desc(-log10(pval)))
# jgene <- jsub$gene[[1]]
# jmark <- "H3K4me1"
#   


jgenes <- c("Hlf", "Tead1", "Hoxa9", "Meis1", "Hoxb9", "Hoxb3", "Hoxd4", "Adgrg1", "S100a7a", "S100a8", "Ltf", "Hdac4", "Hbb-bs", "Hbb-y", "Sox6", "Ebf1", "Irf4", "Pax5", "Cd180", "Cd38", "Iglv3", "Bach2")

clsts.keep <- c("HSPCs", "lymph", "granu", "eryth")
for (jgene in jgenes){
  for (jmark in jmarks){
    jmat <- dat.imputed.lst[[jmark]]
    cnames.keep <- metadat.lst[[jmark]]$cell
    jmat.filt <- jmat[, cnames.keep]
    
    
    jmeta <- metadat.lst[[jmark]] %>%
      dplyr::rename(xvar = cell.var.within.sum.norm, clst = cluster) %>%
      dplyr::select(cell, xvar, clst) %>%
    mutate(xvar.orig = xvar,
           xvar = max(xvar) - xvar)
    jmeta.ordered <- jmeta[match(colnames(jmat.filt), jmeta$cell), ]
    jmeta.ordered$clst <- factor(jmeta.ordered$clst, levels = clsts.keep)
    
    (rowi <- grep(jgene, rownames(jmat.filt), value = TRUE))
    # (rowi <- sample(rownames(jmat.filt), size = 1))
    
    jexprs <- jmat.filt[rowi, ]
    fit.input <- data.frame(exprs = jexprs, xvar = jmeta.ordered$xvar, clst = jmeta.ordered$clst, cell = jmeta.ordered$cell, stringsAsFactors = TRUE) %>%
      rowwise() %>%
      mutate(plate = GetCondFromSamp(samp = cell, mark = jmark))
    
    # addd fit?
    jfit.tmp <- FitClstVar(fit.input$exprs, fit.input)
    exprs.pred <- predict(jfit.tmp)
    fit.input$exprs.pred <- exprs.pred
    
    # # fit null
    # jfit.tmp.null <- FitClstVar.null(fit.input$exprs, fit.input)
    # exprs.pred.null <- predict(jfit.tmp.null)
    # fit.input$exprs.pred.null <- exprs.pred.null
    
    # jname <- jgene
    # jfit <- jfit.tmp
    # jint <- coefficients(jfit)[["(Intercept)"]]  # HSPC intercept
    # # rnames.keep <- which(startsWith(x = names(coefficients(jfit)), prefix = "clst"))
    # rnames.keep <- names(coefficients(jfit))[which(startsWith(x = names(coefficients(jfit)), prefix = "clst"))]
    # jmeans <- jint + coefficients(jfit)[rnames.keep]
    # jslopes <- coefficients(jfit)[which(startsWith(x = names(coefficients(jfit)), prefix = "xvar"))]
    # dat.params <- data.frame(params.mean = c("clstHSC", names(jmeans)), jmean = c(jint, jmeans), params.slope = names(jslopes), jslope = jslopes, stringsAsFactors = FALSE)
    # dat.params$rname <- jname
    # pval.vec.all <- summary(jfit)$coefficients[, "Pr(>|t|)"]
    # pval.dat.keep <- data.frame(params.mean = c("clstHSC", rnames.keep), pval = c(NA, pval.vec.all[rnames.keep]), stringsAsFactors = FALSE)
    # dat.params <- left_join(dat.params, pval.dat.keep)
    # # dat.params$pval <- summary(jfit)$coefficients[, "Pr(>|t|)"]
    
    m.byclst <- ggplot(fit.input, aes(x = xvar, y = exprs, color = clst)) + geom_point(alpha = 0.3)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, rowi)) +  
      geom_line(mapping = aes(x = xvar, y = exprs.pred, color = clst)) + 
      xlab("MNase Overdigestion") + ylab("Imputed ChIC signal")
    
    
    m.byclst.split <- ggplot(fit.input, aes(x = xvar, y = exprs, color = clst)) + geom_point(alpha = 0.3)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, rowi))  + facet_wrap(~clst) + 
      geom_line(mapping = aes(x = xvar, y = exprs.pred), color = "grey55") + 
      xlab("MNase Overdigestion") + ylab("Imputed ChIC signal")
    
    m.bycond <- ggplot(fit.input, aes(x = xvar, y = exprs, color = plate)) + geom_point(alpha = 0.3)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_color_manual(values = cbPalette) + ggtitle(paste(jmark, rowi)) + facet_wrap(~clst) + 
      geom_line(mapping = aes(x = xvar, y = exprs.pred), color = "grey55") + 
      xlab("MNase Overdigestion") + ylab("Imputed ChIC signal")
    
    print(m.byclst)
    # print(m.byclst2)
    print(m.byclst.split)
    print(m.bycond)
  }
}
if (make.plots){
  dev.off()
}


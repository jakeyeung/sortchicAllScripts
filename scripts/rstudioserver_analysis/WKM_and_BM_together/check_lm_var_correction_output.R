# Jake Yeung
# Date of Creation: 2020-06-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/check_lm_var_correction_output.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

FitClstVar <- function(jrow, jmeta){
  fit.input <- data.frame(exprs = jrow, xvar = jmeta$xvar, clst = jmeta$clst, stringsAsFactors = TRUE)
  jfit <- lm(exprs ~ xvar:clst + clst, fit.input)
  return(jfit)
}

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/pdfs_from_WKM_BM_merged"

outf <- file.path(outdir, paste0("plot_varcorrection_and_logfc_examples.", Sys.Date(), ".pdf"))

pdf(file = outf, useDingbats = FALSE)

# Load gene sets ----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

# Load inputs -------------------------------------------------------------

# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits.2020-06-25.RData"
inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/var_slope_estimates/MouseBM_log_lm_fits.2020-06-26.CleanUpEryth.RData"
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


jgenes <- c("Hlf", "Tead1", "Hoxa9", "Meis1", "Hoxb9", "Hoxb3", "Hoxd4", "Adgrg1", "S100a7a", "S100a8", "Ltf", "Hdac4", "Hbb-bs", "Hbb-y", "Sox6", "Ebf1", "Irf4", "Pax5", "Cd180", "Cd38", "Iglv3", "Bach2")

for (jgene in jgenes){
  for (jmark in jmarks){
    jmat <- dat.imputed.lst[[jmark]]
    cnames.keep <- metadat.lst[[jmark]]$cell
    jmat.filt <- jmat[, cnames.keep]
    
    
    jmeta <- metadat.lst[[jmark]] %>%
      # dplyr::rename(xvar = cell.var.within.sum.norm, clst = cluster) %>%
      dplyr::select(cell, xvar, clst)
    # mutate(xvar.orig = xvar,
    #        xvar = max(xvar) - xvar)
    jmeta.ordered <- jmeta[match(colnames(jmat.filt), jmeta$cell), ]
    
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
    print(m.byclst.split)
    print(m.bycond)
  }
}

dev.off()



# Jake Yeung
# Date of Creation: 2020-12-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/10-DE_downstream_all_marks.K9me3_bins.R
# Do DE downstream 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(DescTools)

keeptop <- 150
low.in.k9 <- FALSE
# low.in.k9 <- TRUE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Get DE outputs ----------------------------------------------------------

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
  return(params.dat.all)
})

pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- jfits.lst.lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(jfits.lst.lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < 1e-10)


# Get k9 bins and plot  ---------------------------------------------------


pval.k9.sub <- subset(pvals.lst$H3K9me3, pval < 1e-10) %>%
  arrange(desc(pval))

k9.bins.names <- unique(pval.k9.sub$bin)
ctypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", ctypes.keep, ".Estimate", sep = "")

params.dat.wide <- data.table::dcast(subset(params.lst$H3K9me3, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate") %>%
  rowwise() %>%
  mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
         ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
         ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
         ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
         Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
         Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
         Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
         HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))

params.dat.wide.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(params.lst[[jmark]], bin %in% k9.bins.names & param %in% params.keep) %>%
                   group_by(bin) %>% filter(max(abs(estimate)) < 5)
  jdat <- data.table::dcast(jsub, 
                            formula = bin ~ param, value.var = "estimate") %>%
    rowwise() %>%
    mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
           ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
           ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
           ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
           ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
           ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
           Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
           Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
           Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
           HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))
  # keep only effect cnames 
  cnames.keep.i <- grep("effect$", colnames(jdat))
  cnames.new <- paste(colnames(jdat)[cnames.keep.i], jmark, sep = "_")
  colnames(jdat)[cnames.keep.i] <- cnames.new
  cnames.keep.bin.i <- grep("bin", colnames(jdat))
  cnames.keep.merged.i <- c(cnames.keep.bin.i, cnames.keep.i)
  jdat.filt <- jdat[, cnames.keep.merged.i]
  return(jdat.filt)
})




if (low.in.k9){
  jsort.hspcs <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(HSPCs.effect)
    arrange(desc(HSPCs.effect))
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
  
  jsort.bcell <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(desc(Bcells.effect)) 
    arrange(Bcells.effect)
  jbins.bcell <- jsort.bcell$bin[1:keeptop]
  
  jsort.granu <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(desc(Granulocytes.effect))
    arrange(Granulocytes.effect)
  jbins.granu <- jsort.granu$bin[1:keeptop]
  
  jsort.eryth <- params.dat.wide %>%
    group_by(bin) %>%
    # arrange(descEryths.effect)) 
    arrange(Eryths.effect)
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
} else {
  jsort.hspcs <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(HSPCs.effect)
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]
  
  jsort.bcell <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(desc(Bcells.effect))
  jbins.bcell <- jsort.bcell$bin[1:keeptop]
  
  jsort.granu <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(desc(Granulocytes.effect))
  jbins.granu <- jsort.granu$bin[1:keeptop]
  
  jsort.eryth <- params.dat.wide %>%
    group_by(bin) %>%
    arrange(desc(Eryths.effect))
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
}



bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)

bins.keep.lst <- list("Eryths" = jbins.eryth,
                      "Bcells" = jbins.bcell,
                      "Granulocytes" = jbins.granu,
                      "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

# annotate genes



ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
ctypes.k9me3 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")
dat.metas <- lapply(jmarks, function(jmark){
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/umaps_final/umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt")
  fread(inf)
})

dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})


jmetas.pretty.lst <- lapply(jmarks, function(jmark){
  jmeta <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes)
  } else { 
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes.k9me3)
  }
  jmeta <- jmeta %>% arrange(cluster, jrep)
})

cells.keep.lst <- lapply(jmetas.pretty.lst, function(jdat){
  jdat$cell
})
ctype2col <- hash::hash(jmetas.pretty.lst$H3K4me1$cluster, jmetas.pretty.lst$H3K4me1$colorcode)
names(bins.keep) <- c(rep("Eryths", keeptop), rep("Bcells", keeptop), rep("Granulocytes", keeptop), rep("HSPCs", keeptop))
colsvec <- sapply(names(bins.keep), function(x) AssignHash(x, jhash = ctype2col, null.fill = NA))
bin2col <- hash::hash(bins.keep, colsvec)


# Plot celltype effects label bins ----------------------------------------

# compare ctype effects between two marks (relative to K9me3)
jmerged <- Reduce(f = left_join, x = params.dat.wide.lst[c("H3K4me1", "H3K4me3", "H3K27me3")], init = params.dat.wide.lst$H3K9me3)




# Plot effects: k9me3 on x axis  ------------------------------------------


jmark <- "H3K27me3"
ggplot(params.lst[[jmark]] %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) +
  geom_density() + 
  facet_wrap(~param) + 
  theme_bw() + 
  ggtitle(jmark) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# parms.keep

celltypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", celltypes.keep, ".Estimate", sep = "")

dat.params.wide.lst <- lapply(params.lst, function(jdat){
  jdat.tmp <- jdat %>% 
    filter(param %in% params.keep) %>% 
    group_by(bin) %>%
    filter(max(abs(estimate)) < 5)
  reshape2::dcast(data = jdat.tmp, formula = bin ~ param + mark, value.var = "estimate")
})
dat.params.wide.joined <- Reduce(f = left_join, x = dat.params.wide.lst[c("H3K4me1", "H3K4me3", "H3K27me3")], init = dat.params.wide.lst[["H3K9me3"]])


# jmarkref <- "H3K9me3"
# jmarkref <- "H3K4me1"


# Check Bcell vs Granu estimate  ------------------------------------------



# Load RData  -------------------------------------------------------------

inf.mat.adj <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_H3K4me1_H3K9me3_analysis/batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2020-12-20.H3K27me3rep2rep3reseq.RData"
load(inf.mat.adj, v=T)

mat.adj.lst <- lapply(mat.adj.lst, function(jmat){
  rownames(jmat) <- jmat$rname
  jmat$rname <- NULL
  return(jmat)
})

bins.common <- Reduce(f = intersect, x = lapply(mat.adj.lst, rownames))
bins.keep.common <- bins.keep[bins.keep %in% bins.common]

jmat.lst <- lapply(jmarks, function(jmark){
  jmat <- mat.adj.lst[[jmark]][bins.keep.common, cells.keep.lst[[jmark]]]
  print(dim(jmat))
  jmat <- apply(jmat, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))
  jmat <- t(apply(jmat, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
  return(jmat)
})

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/other_marks"

outpdf1 <- file.path(outdir, paste0("H3K9me3_bins_by_pval.gene_gene_correlations.H3K9me3_vs_other_marks_bins_labeled.LowInK9.", low.in.k9, ".", Sys.Date(), ".bugfix.boxplots.pdf"))

pdf(outpdf1, useDingbats = FALSE)

# check heatmap
for (jmark in jmarks){
  heatmap3::heatmap3(jmat.lst[[jmark]], Rowv = NA, Colv = NA, scale = "row", ColSideColors = jmetas.pretty.lst[[jmark]]$clustercol, RowSideColors = sapply(bins.keep.common, AssignHash, jhash = bin2col),  revC = TRUE, main = paste0(jmark, " 50kb bins"), margins = c(5, 8))
}
 
# make boxplots
for (jmark in jmarks){
  for (jset in bnames){
    m.boxplot <- ggplot(params.lst[[jmark]] %>% 
             filter(bin %in% bins.keep.lst[[jset]]) %>%
             group_by(bin) %>% 
             filter(max(abs(estimate)) < 5), 
           aes(x = forcats::fct_reorder(.f = param, .x = estimate, .fun = median, .desc = TRUE), y = estimate)) + 
      ggtitle(jset, jmark) + 
      ylab("log2FC relative to HSPCs") + 
      geom_boxplot() + 
      geom_point() + 
      theme_bw() + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    print(m.boxplot)
  }
}



# replot but highlight eryth, bcell, granu, and hspc-specific genes

jmarkref <- "H3K9me3"
jmarksfilt <- jmarks[jmarks != jmarkref]
print(jmarksfilt)

for (jmarktmp in jmarksfilt){
  # plot HSPC effect K9me3 vs K4me1
  m <- ggplot(jmerged, aes_string(x = paste0("HSPCs.effect_", jmarkref), y = paste0("HSPCs.effect_", jmarktmp))) + 
    geom_point(alpha = 0.25) +  
    geom_density_2d(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 1.5) + 
    geom_vline(xintercept = 0, linetype = "solid", color = "red", size = 1.5) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(jmerged %>% mutate(in.set = bin %in% bins.keep.lst[["HSPCs"]]), 
              aes_string(x = paste0("HSPCs.effect_", jmarkref), y = paste0("HSPCs.effect_", jmarktmp))) + 
    geom_point(alpha = 0.25, aes(color = in.set)) +  
    geom_density_2d(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 1.5) + 
    geom_vline(xintercept = 0, linetype = "solid", color = "red", size = 1.5) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
} 

for (bname1 in bnames){
  jbins <- bins.keep.lst[[bname1]]
  for (bname2 in bnames){
    print(bname2)
    mlst <- lapply(jmarksfilt, function(jmark){
      jsub <- jmerged %>% 
        mutate(in.set = bin %in% jbins) %>% 
        arrange(in.set)
      m <- ggplot(jsub, 
                  aes_string(x = paste0(bname2, ".effect_", jmarkref), y = paste0(bname2, ".effect_", jmark))) + 
        geom_point(alpha = 0.8, mapping = aes(color = in.set)) + 
        theme_bw() + 
        ggtitle(paste0(bname1, "-spec binslabeled")) + 
        geom_hline(yintercept = 0, linetype = "dotted") + 
        geom_vline(xintercept = 0, linetype = "dotted") + 
        geom_density_2d(alpha = 0.1) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              legend.position = "bottom")
      return(m)
    })
    JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], cols = 3)
  }
}
dev.off()


for (jmarkref in jmarks){
  outpdf <- file.path(outdir, paste0("H3K9me3_bins_by_pval.gene_gene_correlations.reference_", jmarkref, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
  
  mall <- ggplot(params.lst[[jmarkref]] %>% filter(abs(estimate) < 5), aes(x = estimate, fill = param)) +
    geom_density() + 
    facet_wrap(~param) + 
    theme_bw() + 
    ggtitle(jmarkref) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(mall)
  
  
  jmarksfilt <- jmarks[jmarks != jmarkref]
  print(jmarksfilt)
  for (jparam in params.keep){
    mlst <- lapply(jmarksfilt, function(jmark){
      m <- ggplot(dat.params.wide.joined %>% filter(bin %in% k9.bins.names), 
                  aes_string(x = paste0(jparam, "_", jmarkref), y = paste0(jparam, "_", jmark))) + 
        geom_point() + 
        geom_density_2d() + 
        theme_bw() + 
        ggtitle(jmark, jparam) +
        geom_hline(yintercept = 0, linetype = "dotted") + 
        geom_vline(xintercept = 0, linetype = "dotted") + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      return(m)
    })
    JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], cols = 3)
  }
  dev.off()
}


# Write bins to output ----------------------------------------------------






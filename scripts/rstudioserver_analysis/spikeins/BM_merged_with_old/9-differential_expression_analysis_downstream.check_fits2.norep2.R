# Jake Yeung
# Date of Creation: 2020-11-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/8-differential_expression_analysis_downstream.check_fits.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


library(topicmodels)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

jgenes <- c("Pax5", "S100a8", "Hlf", "Tal1", "Hbb-y")
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load fits ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3")

jtype <- "TSS"
for (jmark in jmarks){
  
  # jmark <- "H3K4me1"
  # jmark <- "H3K4me3"
  # jmark <- "H3K27me3"
  outf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/differential_analysis_BM_AllMerged3.gene_examples/DE_downstream_", jtype, ".", jmark, ".with_examples.norep2.pdf")
  
  pdf(file = outf, useDingbats = FALSE)
  
  # inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_TSS_10000.H3K4me1.2020-11-14.newannot.firt_try.RData"
  # load(inf, v=T)
  
  infrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression/integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_TRUE.forPoissonRegression.CountR1only.2020-06-05.smaller.RData"
  assertthat::assert_that(file.exists(infrdata))
  load(infrdata, v=T)
  
  jannots <- names(de.ens.sorted.stringent)
  names(jannots) <- jannots
  
  inf1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_TSS_10000.", jmark, ".2020-11-14.newannot2.norep2.RData")
  inf2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.spikeins/poisson_fit_TSS.dist_10000.", jmark, ".2020-11-14.spikeins.again.newannot2.norep2.RData")
  
  assertthat::assert_that(file.exists(inf1))
  assertthat::assert_that(file.exists(inf2))
  
  out.lst <- lapply(list(total = inf1, spikeins = inf2), FUN = function(inf){
    print(inf)
    load(inf, v=T)
    jrow.names <- names(jfits.lst)
    names(jrow.names) <- jrow.names
    
    # jrow <- grep(jgene, jrow.names, value = TRUE)
    dat.summary <- lapply(jrow.names, function(jrow){
      ctype.effects <- grep("^Cluster", names(jfits.lst[[jrow]]), value = TRUE)
      fits.sub <- unlist(jfits.lst[[jrow]][ctype.effects])
      dat.summary.tmp <- data.frame(param = names(fits.sub), value = fits.sub, rname = jrow, stringsAsFactors = FALSE)
      return(dat.summary.tmp)
    }) %>%
      bind_rows()
    
    dat.summary <- dat.summary %>%
      rowwise() %>%
      mutate(gene = strsplit(rname, "\\.\\.")[[1]][[2]],
             ens = AssignHash(x = gene, jhash = g2e, null.fill = gene))
    return(list(dat.summary = dat.summary, dat.annots.filt.mark = dat.annots.filt.mark, jmat.mark = jmat.mark, ncuts.for.fit.mark = ncuts.for.fit.mark))
  })
  
  
  # Summarize both ---------------------------------------------------------
  
  
  jnames <- names(out.lst)
  names(jnames) <- jnames
  mall.lst <- lapply(jnames, function(jname){
    out <- out.lst[[jname]]
    mall <- ggplot(out$dat.summary, aes(x = value/ log(2), fill = param)) + 
      geom_density(alpha = 0.25) + 
      scale_fill_manual(values = cbPalette) + 
      theme_bw() + 
      facet_wrap(~param) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      ggtitle(jmark, paste0("All Genes. Norm by:", jname)) + 
      coord_cartesian(xlim = c(-10, 10)) + 
      xlab("log2FC relative to HSPCs") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(mall)
  })
  print(mall.lst)
  
  # test some genes
  
  
  for (jgene in jgenes){
    m.byplate.lst <- lapply(jnames, function(jname){
      jylab <- paste0("log2(cuts/", jname, ")")
      dat.summary <- out.lst[[jname]]$dat.summary    
      ncuts.for.fit.mark <- out.lst[[jname]]$ncuts.for.fit.mark
      jmat.mark <- out.lst[[jname]]$jmat.mark
      dat.annots.filt.mark <- out.lst[[jname]]$dat.annots.filt.mark
      jrow <- subset(dat.summary, gene == jgene)$rname[[1]]
      jsub <- subset(dat.summary, gene == jgene)
      
      m.params <- ggplot(jsub, aes(x = param, y = value / log(2))) + 
        geom_point() + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        ylab("log2FC relative to HSPC") + 
        xlab("") + 
        ggtitle(jrow, paste("Norm by:", jname))
      
      # plot log2FCs relative to FC
      dat.counts <- data.frame(cell = colnames(jmat.mark), cuts = jmat.mark[jrow, ], stringsAsFactors = FALSE)
      dat.counts.annot <- left_join(dat.annots.filt.mark, dat.counts, by = "cell") %>%
        left_join(., ncuts.for.fit.mark, by = "cell") %>%
        rowwise() %>%
        mutate(Experi = ClipLast(x = Plate, jsep = "-"))
      
      m <- ggplot(dat.counts.annot, aes(x = cluster, y = log2(cuts / ncuts.total))) + 
        facet_wrap(~Plate) + 
        geom_boxplot() + 
        geom_point() +
        ggtitle(paste(jrow), paste0("Norm by:", jname)) + 
        theme_bw() + 
        ylab(jylab) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      mlst <- list(m.params, m)
      return(mlst)
    })
    print(m.byplate.lst)
  }
  
  dev.off()
  
  
}

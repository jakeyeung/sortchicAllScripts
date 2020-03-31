# Jake Yeung
# Date of Creation: 2020-03-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/multiomics_integration.R
# Multi-omics integration of celltypes 

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

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#fa8c30", "#bd0a42", "#1347d5")
jdists <- c(500L, 1000L, 5000L, 10000L)
jsuffixs <- c(".celltypes_filt", "")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
m2c <- MarkerToCelltype()

make.plots <- TRUE
# jdist <- 1000L
# jsuffix <- jsuffixs[[1]]

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/proms_enhs_genebody_analysis"
dir.create(outdir)

for (jdist in jdists){
  for (jsuffix in jsuffixs){
   
     
    print(Sys.time() - jstart)
    # jdist <- 1000L
    # jsuffix <- ".celltypes_filt"
    outf <- file.path(outdir, paste0("multiomics_summary_ctypes.dist_", jdist, jsuffix, ".", Sys.Date(), ".pdf"))
    
    if (file.exists(outf)){
      print(paste("outf exists, skipping:", outf))
      next
    }
    
    if (make.plots){
      # width=1920&height=989. 72 pixels per inch when you click zoom in Rstudio 
      pdf(outf, useDingbats = FALSE, width = 1920 / 72, height = 989 / 72)
    }
    
    inf.de <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat", jsuffix, ".rds")
    assertthat::assert_that(file.exists(inf.de))
    # inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.celltypes_filt.rds"
    # inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Giladi_et_al_2018/diff_exprs_Giladi_seurat.rds"
        
    # get raw counts ----------------------------------------------------------
    
    indir.mat <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.AnnotatedGeneRegionsWithPromsEnhs")
    count.mat.lst.pe <- lapply(jmarks, function(jmark){
      print(jmark)
      inf.mat <- file.path(indir.mat, paste0(jmark, ".countTableTSS.mapq_40.TSS_", jdist, ".blfiltered.csv"))
      jmat <- ReadMatTSSFormat(inf.mat, as.sparse = TRUE, add.coord = TRUE)
      cols.i <- order(colnames(jmat))
      jmat <- jmat[, cols.i]
    })
    
    
    # add gene bodies
    indir.mat.gb <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.AnnotatedGeneRegionsWithGeneBodies")
    assertthat::assert_that(dir.exists(indir.mat.gb))
    
    count.mat.lst.gb <- lapply(jmarks, function(jmark){
      print(jmark)
      inf.mat <- file.path(indir.mat.gb, paste0(jmark, ".countTableTSS.mapq_40.gene_body.blfiltered.csv"))
      assertthat::assert_that(file.exists(inf.mat))
      jmat <- ReadMatTSSFormat(inf.mat, as.sparse = TRUE, add.coord = TRUE)
      cols.i <- order(colnames(jmat))
      jmat <- jmat[, cols.i]
    })
    
    # merge with promote and enhancer matrix??
    
    count.mat.lst <- lapply(jmarks, function(jmark){
      jmat1 <- count.mat.lst.pe[[jmark]]
      jmat2 <- count.mat.lst.gb[[jmark]]
      assertthat::assert_that(identical(colnames(jmat1), colnames(jmat2)))
      # simple rbind
      jmat.merge <- rbind(jmat1, jmat2)
    })
    
    # check rownames
    rnames.all <- lapply(count.mat.lst, function(x) rownames(x))
    lapply(rnames.all, length)
    
    rnames.common <- Reduce(intersect, rnames.all)
    
    print(length(rnames.common))
    
    count.mat.lst.filt <- lapply(count.mat.lst, function(x){
      x[rnames.common, ]
    })
    
    lapply(count.mat.lst.filt, dim)
    
    # label each row as either promoter, enhancer, or genebody
    geneannot <- hash()
    # get promoter regions
    rnames.pe <- lapply(count.mat.lst.pe, function(x){
      rownames(x)
    }) %>%
      unlist() %>%
      unique()
    rnames.gb <- lapply(count.mat.lst.gb, function(x){
      rownames(x)
    }) %>%
      unlist() %>%
      unique()
    
    rnames.enh <- rnames.pe[grepl("enhancer", rnames.pe)]
    rnames.prom <- rnames.pe[!grepl("enhancer", rnames.pe)]
    
    rnames.enh.lst <- rep("enhancer", length(rnames.enh)); names(rnames.enh.lst) <- rnames.enh
    rnames.prom.lst <- rep("promoter", length(rnames.prom)); names(rnames.prom.lst) <- rnames.prom
    rnames.gb.lst <- rep("genebody", length(rnames.gb)); names(rnames.gb.lst) <- rnames.gb
    
    # # hash it for speed
    # biotype.hash <- hash::hash()
    # for (jlst in list(rnames.enh.lst, rnames.prom.lst, rnames.gb.lst)){
    #   for (jname in names(jlst)){
    #     biotype.hash[[jname]] <- jlst[[jname]]
    #   }
    # }
    biotype.dat <- data.frame(region_coord_full = c(names(rnames.enh.lst), names(rnames.prom.lst), names(rnames.gb.lst)),
                              biotype = c(rnames.enh.lst, rnames.prom.lst, rnames.gb.lst), 
                              stringsAsFactors = FALSE)
    # biotype.hash <- hash(biotype.dat$region_coord_full, biotype.dat$biotype)
    
    
    # Celltype sepcific genes -------------------------------------------------
    
    pval.min <- 0.01
    fc.min <- 0
    dat.de <- readRDS(inf.de)
    neutro.genes <- subset(dat.de, p_val_adj < pval.min & cluster == "Ltf" & avg_logFC > fc.min)$gene
    
    # label genes as cluster specific
    
    jclsts <- as.character(unique(dat.de$cluster))
    names(jclsts) <- jclsts
    
    de.genes.lst <- lapply(jclsts, function(jclst){
      subset(dat.de, p_val_adj < pval.min & cluster == jclst & avg_logFC > fc.min)$gene
    })
    
    subset(dat.de, gene == "S100a8")
    
    
    gene2ctype <- hash::hash()
    for (jclst in jclsts){
      for (jgene in de.genes.lst[[jclst]]){
        gene2ctype[[jgene]] <- c(gene2ctype[[jgene]], jclst)
      }
    }
    
    # Load annots  ------------------------------------------------------------
    
    indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")
    dat.annots.all <- lapply(jmarks, function(jmark){
      inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
      load(inf.annots, v=T)
      return(dat.umap.glm.fillNAs)
    })
    
    # this is slow... make this faster
    outfcellsizes <- file.path("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs", paste0("glmpca_cellsizes_BM_AllMerged.", Sys.Date(), ".rds"))
    if (file.exists(outfcellsizes)){
      dat.totalcuts <- readRDS(outfcellsizes)
    } else {
      indir.totalcuts <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs")
      dat.totalcuts <- lapply(jmarks, function(jmark){
        inf.glmpca <- file.path(indir.totalcuts, paste0("PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData"))
        load(inf.glmpca, v=T)
        return(glm.inits$size.factor)
      })
      saveRDS(dat.totalcuts, file = outfcellsizes)
    }
    
    
    # Define each dot as either TSS or enhancer -------------------------------
    
    coords.vec <- sapply(rnames.common, function(x) strsplit(x, ";")[[1]][[1]])
    # what to do with TSS of genes with exact same TSS? remove it? yeah... they should get same label
    # e.g. chrX:135732733-135752733;Armcx5+ chrX:135732733-135752733;Gprasp1
    coords.vec <- coords.vec[!duplicated(coords.vec)]
    jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/first_transcript_tss/gene_tss_winsize.50000.first_transcript.bed"
    # annotate coord to gene
    jchromos.keep <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
    bins.annot <- AnnotateCoordsFromList(coords.vec = coords.vec, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos.keep)
    
    # Make pseudobulks  -------------------------------------------------------
    
    cnames.keep.lst.all <- lapply(jmarks, function(jmark){
      jsplit <- split(dat.annots.all[[jmark]], dat.annots.all[[jmark]]$cluster)
      cnames.keep <- lapply(jsplit, function(x) x$cell)
    })
    # merge clusters
    
    # clsts.merge <- list(H3K4me1 = hash::hash(c("Bcells-Cd47_topic29", "Bcells-Cd83_topic10"), "BcellsMerged"),
    #                     H3K4me3 = hash::hash(c("Eryth-Cdk6_topic9", "Eryth-Gfi1b_topic7", "Eryth-Sox6_topic16"), "ErythMerged"),
    #                     H3K27me3 = hash::hash(c("Eryth-Gfi1-_topic17", "Eryth-Slc7a6-_topic1", "Eryth-Sox6-_topic6"), "ErythMerged"))
    
    clsts.merge <- list(H3K4me1 = list("BcellsMerged" = c("Bcells-Cd47_topic29", "Bcells-Cd83_topic10"),
                                       "InnateLymphMerged" = c("ILC-PrkchPlus_topic18", "ILC-RoraPlus_topic11")),
                        H3K4me3 = list("ErythMerged" = c("Eryth-Cdk6_topic9", "Eryth-Gfi1b_topic7", "Eryth-Sox6_topic16")),
                        H3K27me3 = list("ErythMerged" = c("Eryth-Gfi1-_topic17", "Eryth-Slc7a6-_topic1", "Eryth-Sox6-_topic6")))
    
    cnames.keep.lst.all <- lapply(jmarks, function(jmark){
      for (newname in names(clsts.merge[[jmark]])){
        assertthat::assert_that(length(newname) == 1)
        oldnames <- clsts.merge[[jmark]][[newname]]
        cnames.merge <- unique(unlist(cnames.keep.lst.all[[jmark]][oldnames]))
        cnames.keep.lst.all[[jmark]][[newname]] <- cnames.merge
        # remove old names
        for (oldname in oldnames){
          cnames.keep.lst.all[[jmark]][[oldname]] <- NULL
        }
      }
      return(cnames.keep.lst.all[[jmark]])
    })
    # lapply(cnames.keep.lst.all, function(x) names(x))
    
    
    # Get pseudobulks ---------------------------------------------------------
    
    
    count.pseudos <- lapply(jmarks, function(jmark){
      exprs.lst <- SumAcrossClusters(count.mat.lst.filt[[jmark]], cnames.keep.lst.all[[jmark]])
      exprs.mat <- as.data.frame(do.call(cbind, exprs.lst))
      cols.i <- order(colnames(exprs.mat))
      exprs.mat[, cols.i]
    })
    
    
    # Create renaming list ----------------------------------------------------
    
    
    
    lapply(count.pseudos, function(x) colnames(x))
    
    
    # filter out relevant celltypes and create "2D" matrix?
    h3k4me1.cnames <- list("Neutrophils_topic23" = "Neutro", 
                           "HSCs-Hlf_topic7" = "HSCs",
                           "Eryth_topic27" = "Eryth",
                           "InnateLymphMerged" = "NKcells",
                           # "Bcells-Cd47_topic29" = "Bcells",
                           "BcellsMerged" = "Bcells")
    
    h3k4me3.cnames <- list("Bcells_topic13" = "Bcells", 
                           # "Eryth-Sox6_topic16" = "Eryth",
                           "ErythMerged" = "Eryth", 
                           "HSCs-Hlf_topic26" = "HSCs",
                           "Neutrophils_topic2" = "Neutro",
                           "InnateLymph_topic27" = "NKcells")
    
    h3k27me3.cnames <- list("Bcells_topic16" = "Bcells", 
                            # "Eryth-Sox6-_topic6" = "Eryth",
                            "ErythMerged" = "Eryth",
                            "HSCs-Tead1-_topic9" = "HSCs",
                            "Neutrophils_topic22" = "Neutro",
                            "InnateLymph_topic27" = "NKcells")
    print(jmarks)
    jmarks.cnames <- list(H3K4me1 = h3k4me1.cnames,
                          H3K4me3 = h3k4me3.cnames,
                          H3K27me3 = h3k27me3.cnames)
    
    
    
    # Rename column names, collapse names with same names ---------------------
    
    count.mat.pbulk.all <- lapply(jmarks, function(jmark){
      cnames <- jmarks.cnames[[jmark]]
      count.mark <- lapply(names(cnames), function(cname){
        x <- count.pseudos[[jmark]][[cname]]
        assertthat::assert_that(!is.null(x))
        names(x) <- rnames.common
        return(x)
      })
      names(count.mark) <- names(cnames)
      # rename
      names(count.mark) <- sapply(names(cnames), function(x) cnames[[x]])
      count.mat <- do.call(cbind, count.mark)
      # order order column names alphabetically, take common rnames
      cols.i <- order(colnames(count.mat))
      count.mat <- count.mat[rnames.common, cols.i]
    })
    
    
    
    # get normalization constant by summing up across cells, also track how many cells per cluster ----
    
    totalcuts.pseudos <- lapply(jmarks, function(jmark){
      exprs.lst <- lapply(cnames.keep.lst.all[[jmark]], function(cnames){
        ikeep <- which(names(dat.totalcuts[[jmark]]) %in% cnames)
        x <- dat.totalcuts[[jmark]][which(names(dat.totalcuts[[jmark]]) %in% cnames)]
        totalcuts <- sum(x)
      })
      exprs.mat <- as.data.frame(do.call(cbind, exprs.lst))
      cnames.i <- order(colnames(exprs.mat))
      exprs.mat <- exprs.mat[, cnames.i]
      return(exprs.mat)
    })
    
    ncells.pseudos <- lapply(jmarks, function(jmark){
      exprs.lst <- lapply(cnames.keep.lst.all[[jmark]], function(cnames){
        return(length(cnames))
      })
      exprs.mat <- as.data.frame(do.call(cbind, exprs.lst))
    })
    
    totalcuts.dat <- lapply(jmarks, function(jmark){
      x <- totalcuts.pseudos[[jmark]]
      dat <- data.frame(pseudobulk = names(x), totalcuts = unlist(x), stringsAsFactors = FALSE)
      dat$mark <- jmark
      return(dat)
    })
    
    ncells.dat <- lapply(jmarks, function(jmark){
      x <- ncells.pseudos[[jmark]]
      dat <- data.frame(pseudobulk = names(x), ncells = unlist(x), stringsAsFactors = FALSE)
      dat$mark <- jmark
      return(dat)
    })
    
    mlst <- lapply(totalcuts.dat, function(dat){
      m <- ggplot(dat, aes(x = pseudobulk, y = totalcuts)) + geom_col() + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        facet_wrap(~mark) 
      return(m)
    })
    JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], cols = 3)
    
    mlst <- lapply(ncells.dat, function(dat){
      m <- ggplot(dat, aes(x = pseudobulk, y = ncells)) + geom_col() + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        facet_wrap(~mark) 
      return(m)
    })
    JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], cols = 3)
    
    ncellscounts.dat <- lapply(jmarks, function(jmark){
      left_join(ncells.dat[[jmark]], totalcuts.dat[[jmark]])
    })
    
    mlst <- lapply(ncellscounts.dat, function(dat){
      m <- ggplot(dat, aes(x = pseudobulk, y = totalcuts / ncells)) + geom_col() + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        facet_wrap(~mark) 
      return(m)
    })
    JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], cols = 3)
    
    
    # Normalize by library size  ----------------------------------------------
    
    totalcuts.pseudos.renamed <- lapply(jmarks, function(jmark){
      x <- totalcuts.pseudos[[jmark]]
      jnames <- jmarks.cnames[[jmark]]
      xfilt <- x[which(names(x) %in% names(jnames))]
      names(xfilt) <- sapply(names(xfilt), function(y) jnames[[y]])
      return(xfilt[order(names(xfilt))])
    })
    
    size.factors.lst <- lapply(jmarks, function(jmark){
      tc <- totalcuts.pseudos.renamed[[jmark]]  # column order shold match colnames in count.mat
      assertthat::assert_that(identical(names(tc), colnames(count.mat.pbulk.all[[jmark]])))
      return(tc)
    })
    
    # normalize
    count.mat.pbulk.all.norm.linear <- lapply(jmarks, function(jmark){
      matnorm <- sweep(count.mat.pbulk.all[[jmark]], MARGIN = 2, STATS = unlist(size.factors.lst[[jmark]]), FUN = "/")
    })
    
    count.mat.pbulk.all.norm.linear.allpseudos <- lapply(jmarks, function(jmark){
      matnorm <- sweep(count.pseudos[[jmark]], MARGIN = 2, STATS = unlist(totalcuts.pseudos[[jmark]]), FUN = "/")
    })
    
    # how much proms and enhs cover total counts
    lapply(count.mat.pbulk.all.norm.linear, function(mat) colSums(mat))
    
    # split by promoter and enhancers then plot for now 
    dat.frac.promsenhs <- lapply(jmarks, function(jmark){
      jmat <- count.mat.pbulk.all.norm.linear[[jmark]]
      rows.enhs <- which(rownames(jmat) %in% subset(biotype.dat, biotype == "enhancer")$region_coord_full)
      rows.proms <- which(rownames(jmat) %in% subset(biotype.dat, biotype == "promoter")$region_coord_full)
      rows.gb <- which(rownames(jmat) %in% subset(biotype.dat, biotype == "genebody")$region_coord_full)
      print(jmark) 
      sapply(list(rows.enhs, rows.proms, rows.gb), function(x) print(length(x)))
      frac.enhs <- colSums(jmat[rows.enhs, ])
      frac.proms <- colSums(jmat[rows.proms, ])
      frac.gb <- colSums(jmat[rows.gb, ])
      
      
      assertthat::assert_that(identical(names(frac.enhs), names(frac.proms)))
      all(sapply(list(names(frac.enhs), names(frac.proms)), FUN = identical, names(frac.gb)))
      frac.dat.wide <- data.frame(pseudobulks = names(frac.enhs), frac.enhs = frac.enhs, frac.proms = frac.proms, frac.gb = frac.gb, stringsAsFactors = FALSE)
      frac.dat <- reshape2::melt(frac.dat.wide, id.vars = "pseudobulks", variable.name = c("biotype"), value.name = "FractionCuts")
    })
    
    m.fracs <- lapply(jmarks, function(jmark){
      m <- ggplot(dat.frac.promsenhs[[jmark]], aes(x = pseudobulks, y = FractionCuts, group = biotype, fill = biotype)) + 
        geom_col(position = "dodge") + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        ggtitle(jmark)
      print(m)
    })
    
    # check all pseduboulks
    dat.frac.promsenhs.allpseudos <- lapply(jmarks, function(jmark){
      jmat.unnorm <- count.pseudos[[jmark]]
      sf <- unlist(totalcuts.pseudos[[jmark]])
      assertthat::assert_that(identical(colnames(jmat.unnorm), names(sf)))
      jmat <- sweep(jmat.unnorm, MARGIN = 2, STATS = sf, FUN = "/")
      
      rows.enhs <- which(rownames(jmat) %in% subset(biotype.dat, biotype == "enhancer")$region_coord_full)
      rows.proms <- which(rownames(jmat) %in% subset(biotype.dat, biotype == "promoter")$region_coord_full)
      rows.gb <- which(rownames(jmat) %in% subset(biotype.dat, biotype == "genebody")$region_coord_full)
      print(jmark) 
      sapply(list(rows.enhs, rows.proms, rows.gb), function(x) print(length(x)))
      frac.enhs <- colSums(jmat[rows.enhs, ])
      frac.proms <- colSums(jmat[rows.proms, ])
      frac.gb <- colSums(jmat[rows.gb, ])
      
      frac.dat.wide <- data.frame(pseudobulks = names(frac.enhs), frac.enhs = frac.enhs, frac.proms = frac.proms, frac.gb = frac.gb, stringsAsFactors = FALSE)
      frac.dat <- reshape2::melt(frac.dat.wide, id.vars = "pseudobulks", variable.name = c("biotype"), value.name = "FractionCuts")
    })
    
    m.fracs.allpseudos <- lapply(jmarks, function(jmark){
      print(jmark)
      m <- ggplot(dat.frac.promsenhs.allpseudos[[jmark]], aes(x = pseudobulks, y = FractionCuts, group = biotype, fill = biotype)) + 
        geom_col(position = "dodge") + 
        theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        ggtitle(jmark)
      print(m)
    })
    
    # Do log normalizatoin and continue?  -------------------------------------
    
    m2 <- lapply(jmarks, function(jmark){
      m <- ggplot(totalcuts.dat[[jmark]], aes(x = pseudobulk, y = totalcuts)) + 
        geom_col() + theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      return(m)
    })
    m2
    
    count.mat.pbulk.all.norm <- lapply(count.mat.pbulk.all.norm.linear, function(x){
      xlog <- log2(x * 10^6 + 1)
      # remove rows that are all zero
      rows.keep <- which(rowSums(xlog) > 0)
      xlog.filt <- xlog[rows.keep, ]
      return(xlog.filt)
    })
    
    count.mat.pbulk.all.norm.allpseudos <- lapply(count.mat.pbulk.all.norm.linear.allpseudos, function(x){
      xlog <- log2(x * 10^6 + 1)
      # remove rows that are all zero
      rows.keep <- which(rowSums(xlog) > 0)
      xlog.filt <- xlog[rows.keep, ]
      return(xlog.filt)
    })
    
    # reset rnames common again?? maybe only if necessary
    rnames2 <- lapply(count.mat.pbulk.all.norm, function(x) rownames(x))
    rnames.common2 <- Reduce(intersect, rnames2)
    
    lapply(jmarks, function(jmark){
      boxplot(data.frame(count.mat.pbulk.all.norm[[jmark]]), las = 2, main = jmark)
    })
    
    lapply(jmarks, function(jmark){
      boxplot(data.frame(count.mat.pbulk.all.norm.allpseudos[[jmark]]), las = 2, main = jmark)
    })
    
    
    # Annotate bins -----------------------------------------------------------
    
    count.mat.pbulk.all.norm.wide <- lapply(jmarks, function(jmark){
      cmat <- count.mat.pbulk.all.norm[[jmark]][rnames.common2, ]
      colnames(cmat) <- paste(jmark, colnames(cmat), sep = "_")
      return(cmat)
    }) 
    count.mat.pbulk.all.norm.wide <- do.call(cbind, count.mat.pbulk.all.norm.wide)
    
    counts.pbulk.long <- data.table::melt(data.frame(region_coord_full = rownames(count.mat.pbulk.all.norm.wide), 
                                                     count.mat.pbulk.all.norm.wide, stringsAsFactors = FALSE), 
                                          id.vars = "region_coord_full", variable.name = "pseudobulk", value.name = "exprs") %>%
      ungroup() %>% 
      mutate(region_coord = sapply(region_coord_full, function(x) strsplit(x, split = ";")[[1]][[1]]),
             mark = sapply(as.character(pseudobulk), function(x) strsplit(x, "_")[[1]][[1]]),
             ctype = sapply(as.character(pseudobulk), function(x) strsplit(x, "_")[[1]][[2]])) %>%
      left_join(., biotype.dat) %>%
      left_join(., bins.annot$out2.df.closest) %>%
      mutate(abs.dist.to.tss = abs(dist.to.tss)) %>%
      filter(!is.na(gene))
    
    counts.pbulk.long.lst <- lapply(jmarks, function(jmark){
      jmat <- count.mat.pbulk.all.norm[[jmark]]  # need to make Bcells -> H3K4me3_Bcells for downstream??
      jlong <- data.table::melt(data.frame(region_coord_full = rownames(jmat), 
                                           jmat, stringsAsFactors = FALSE), 
                                id.vars = "region_coord_full", variable.name = "pseudobulk", value.name = "exprs") %>%
        ungroup() %>%
        mutate(region_coord = sapply(region_coord_full, function(x) strsplit(x, split = ";")[[1]][[1]]),
               mark = jmark,
               ctype = as.character(pseudobulk)) %>%
        left_join(., biotype.dat) %>%
        left_join(., bins.annot$out2.df.closest) %>%
        mutate(abs.dist.to.tss = abs(dist.to.tss)) %>%
        filter(!is.na(gene))
      return(jlong)
    })
    
    counts.pbulk.lst.allpseudos <- lapply(jmarks, function(jmark){
      jmat <- count.mat.pbulk.all.norm.allpseudos[[jmark]]  # need to make Bcells -> H3K4me3_Bcells for downstream??
      jlong <- data.table::melt(data.frame(region_coord_full = rownames(jmat), 
                                           jmat, stringsAsFactors = FALSE), 
                                id.vars = "region_coord_full", variable.name = "pseudobulk", value.name = "exprs") %>%
        ungroup() %>%
        mutate(region_coord = sapply(region_coord_full, function(x) strsplit(x, split = ";")[[1]][[1]]),
               mark = jmark,
               ctype = as.character(pseudobulk)) %>%
        left_join(., biotype.dat) %>%
        left_join(., bins.annot$out2.df.closest) %>%
        mutate(abs.dist.to.tss = abs(dist.to.tss)) %>%
        filter(!is.na(gene))
      return(jlong)
    })
    counts.pbulk.long.allpseudos <- counts.pbulk.lst.allpseudos %>% bind_rows()
    
    
    counts.pbulk.long.allpseudos$de.ctype <- sapply(counts.pbulk.long.allpseudos$gene, AssignHash, gene2ctype, null.fill = NA)
    counts.pbulk.long.allpseudos <- counts.pbulk.long.allpseudos %>%
      group_by(region_coord_full) %>%
      mutate(de.ctype.choose = sample(de.ctype[[1]], size = 1))
    
    counts.pbulk.lst.allpseudos <- lapply(counts.pbulk.lst.allpseudos, function(jdat){
      jdat$de.ctype <- sapply(jdat$gene, AssignHash, gene2ctype, null.fill = NA)
      jdat <- jdat %>%
        group_by(region_coord_full) %>%
        mutate(de.ctype.choose = sample(de.ctype[[1]], size = 1))
    })
    
    counts.pbulk.long$de.ctype <- sapply(counts.pbulk.long$gene, AssignHash, gene2ctype, null.fill = NA)
    counts.pbulk.long <- counts.pbulk.long %>%
      group_by(region_coord_full) %>%
      mutate(de.ctype.choose = sample(de.ctype[[1]], size = 1))
    
    counts.pbulk.long.lst <- lapply(counts.pbulk.long.lst, function(jdat){
      jdat$de.ctype <- sapply(jdat$gene, AssignHash, gene2ctype, null.fill = NA)
      jdat <- jdat %>%
        group_by(region_coord_full) %>%
        mutate(de.ctype.choose = sample(de.ctype[[1]], size = 1))
    })
    
    
    
    # PCA  --------------------------------------------------------------------
    
    # filter for only DE genes?
    # de.genes.lst
    
    de.genes <- unlist(de.genes.lst)
    
    coords.keep <- subset(bins.annot$out2.df.closest, gene %in% de.genes)$region_coord
    
    jfilt <- unique(subset(counts.pbulk.long, region_coord %in% coords.keep)$region_coord_full)
    
    # should be mark by mark....
    pca.out.lst <- lapply(count.mat.pbulk.all.norm.allpseudos, function(jmat){
      #  jmat <- jmat[jfilt, ]
      jmat <- jmat[jfilt, ]
      print(dim(jmat))
      pca.out <- prcomp(t(jmat), center = TRUE, scale. = TRUE)
      dat.pca <- data.frame(sample = rownames(pca.out$x), pca.out$x, stringsAsFactors = FALSE)
      dat.loadings <- data.frame(bin = rownames(pca.out$rotation), pca.out$rotation, stringsAsFactors = FALSE)
      # mutate(bin.labelPC1 <- ifelse(PC1 > quantile(0.99)))
      dat.loadings$bin.labelPC1 <- dat.loadings
      return(list(dat.pca = dat.pca, dat.loadings = dat.loadings))
    })
    
    # plot
    m.pca <- lapply(jmarks, function(jmark){
      out.pca <- pca.out.lst[[jmark]]
      m <- ggplot(out.pca$dat.pca, aes(x = PC1, y = PC2, label = sample)) + 
        geom_point() + geom_text_repel() + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        ggtitle(jmark)
    })
    print(m.pca)
    
    
    # Can we see why the HSCs are high in H3K4me3? -------------------------------
    
    # m.loadings <- lapply(jmarks, function(jmark){
    #   out.pca <- pca.out.lst[[jmark]]
    #   # only label top 10 bins for PC1
    #   m <- ggplot(out.pca$dat.loadings, aes(x = PC1, y = PC2)) + 
    #     geom_point() + 
    #     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    #     ggtitle(jmark)
    # })
    # print(m.loadings)
    
    
    # Confirm celltype specfic signal?  ---------------------------------------
    
    
    
    # counts.pbulk.long$de.ctype.choose <- sapply(counts.pbulk.long$de.ctype, function(x) sample(x = x, size = 1))
    
    subset(counts.pbulk.long, !is.na(ctype))
    dim(subset(counts.pbulk.long, !is.na(ctype)))
    subset(counts.pbulk.long, gene == "S100a8")
    
    # eg chr1:34449762-34469762;Ptpn18
    # subset(counts.pbulk.long, region_coord_full == "chr1:34449762-34469762;Ptpn18+")$de.ctype.choose
    
    # Check celltype specificity  ---------------------------------------------
    
    # on H3K4me3
    # jmark <- "H3K27me3"
    jmark <- "H3K4me1"
    jmark <- "H3K4me3"
    
    # mshapes <- lapply(jmarks, function(jmark){
    #   jsub <- counts.pbulk.long.lst[[jmark]]
    #   m <- ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
    #     facet_grid(ctype ~ de.ctype.choose) + 
    #     # facet_wrap(~biotype, ncol = 1) + 
    #     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    #     scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
    #   return(m)
    # })
    # print(mshapes)
    
    mshapes3 <- lapply(jmarks, function(jmark){
      jsub <- counts.pbulk.long.lst[[jmark]] %>%
        group_by(region_coord_full) %>%
        mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
               logFC = scale(exprs, center = TRUE, scale = FALSE))
      m.exprs <- ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
        facet_grid(ctype ~ de.ctype.choose) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
      m.zscore <- ggplot(jsub, aes(x = zscore, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
        facet_grid(ctype ~ de.ctype.choose) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
      m.logFC <- ggplot(jsub, aes(x = logFC, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
        facet_grid(ctype ~ de.ctype.choose) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
      return(list(m.exprs, m.zscore, m.logFC))
    })
    print(mshapes3)
    
    
    mshapes3.all <- lapply(jmarks, function(jmark){
      jsub <- counts.pbulk.lst.allpseudos[[jmark]] %>%
        group_by(region_coord_full) %>%
        mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
               logFC = scale(exprs, center = TRUE, scale = FALSE))
      m.exprs <- ggplot(jsub, aes(x = exprs, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
        facet_grid(ctype ~ de.ctype.choose) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
      m.zscore <- ggplot(jsub, aes(x = zscore, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
        facet_grid(ctype ~ de.ctype.choose) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
      m.logFC <- ggplot(jsub, aes(x = logFC, fill = biotype, group = biotype)) + geom_density(alpha = 0.3) + 
        facet_grid(ctype ~ de.ctype.choose) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_fill_manual(values = cbPalette) + ggtitle(paste0(jmark, " Bin dist:", jdist))
      return(list(m.exprs, m.zscore, m.logFC))
    })
    print(mshapes3.all)
    
    
    # Plot H3K4me1 vs H3K27me3  -----------------------------------------------
    
    # the random de.ctype.choose is problematic here!
    jctypes <- sort(unique(counts.pbulk.long$ctype))
    for (jctype in jctypes){
      jsub <- subset(counts.pbulk.long, ctype == jctype)
      jsub.wide <- reshape2::dcast(subset(jsub, select = c(region_coord_full, exprs, mark, de.ctype.choose, biotype, gene, de.ctype.choose)), 
                                   formula = "region_coord_full + gene + biotype + de.ctype.choose ~ mark", value.var = "exprs")
      m <- ggplot(jsub.wide %>% arrange(de.ctype.choose) %>% filter(!is.na(de.ctype.choose)), 
             aes(x = H3K4me3, y = H3K27me3, color = de.ctype.choose)) + geom_point() + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jctype) + 
        facet_wrap(~de.ctype.choose)
      print(m)
    }
    # jctype <- "NKcells"
    
    # Compare one ctype with another ------------------------------------------
    
    jsub.wide.lst <- lapply(jmarks, function(jmark){
      jsub <- counts.pbulk.long.lst[[jmark]] %>%
        group_by(region_coord_full)
      jsub.wide <- reshape2::dcast(subset(jsub, select = c(region_coord_full, de.ctype.choose, biotype, gene, de.ctype.choose, ctype, exprs)), 
                                   formula = "region_coord_full  + biotype + gene + de.ctype.choose ~ ctype", value.var = "exprs")
    })
    
    ctypes <- sort(unique(counts.pbulk.long$ctype))
    m.ctypevctype <- lapply(jmarks, function(jmark){
      print(jmark)
      for (ctype1 in ctypes){
        for (ctype2 in ctypes){
          if (ctype1 == ctype2){
            next
          }
          print(paste(ctype1, ctype2))
          m <- ggplot(jsub.wide.lst[[jmark]] %>% arrange(de.ctype.choose) %>% filter(!is.na(de.ctype.choose)), 
                      aes_string(x = ctype1, y = ctype2, color = "biotype")) + geom_point() + 
            theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(jmark) + 
            facet_wrap(~de.ctype.choose)
          print(m)
        }
      }
    })
    
    
    # Compare promoter and enhancer signal in K4me1 vs K4me3  -----------------
    
    stripsize <- 7
    
    jsub <- counts.pbulk.long %>%
      group_by(region_coord_full, mark) %>%
      mutate(zscore = scale(exprs, center = TRUE, scale = TRUE),
             logFC = scale(exprs, center = TRUE, scale = FALSE))
    
    # show H3k4me1, H3K4me3, H3K27me3 levels for each gene? 
    m <- ggplot(jsub, aes(x = ctype, y = exprs, fill = mark)) + geom_boxplot() + 
      theme_bw() + 
      facet_grid(biotype ~ de.ctype.choose) + ggtitle(paste0("Bin size:", jdist)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            strip.text = element_text(size = stripsize, margin = margin()),
            legend.position = "bottom") + 
      xlab("")
    print(m)
    
    m <- ggplot(jsub, aes(x = ctype, y = zscore, fill = mark)) + geom_boxplot() + 
      theme_bw() + 
      facet_grid(biotype ~ de.ctype.choose) + ggtitle(paste0("Bin size:", jdist)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            strip.text = element_text(size = stripsize, margin = margin()),
            legend.position = "bottom") + 
      xlab("")
    print(m)
    
    m <- ggplot(jsub, aes(x = ctype, y = logFC, fill = mark)) + geom_boxplot() + 
      theme_bw() + 
      facet_grid(biotype ~ de.ctype.choose) + ggtitle(paste0("Bin size:", jdist)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            strip.text = element_text(size = stripsize, margin = margin()),
            legend.position = "bottom") + 
      xlab("")
    print(m)
    
    if (make.plots){
      dev.off()
      print(paste("Done writing:", outf))
    }
        
  }
}


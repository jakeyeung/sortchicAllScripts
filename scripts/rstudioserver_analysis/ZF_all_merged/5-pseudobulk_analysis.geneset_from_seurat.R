# Jake Yeung
# Date of Creation: 2020-04-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/5-pseudobulk_analysis.R
# Get pseudobulks from TSS signal???

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(DESeq2)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

EnsemblGene2Gene.ZF <- function(gene.list, return.original=TRUE, species = "drerio"){
  library("biomaRt")
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = paste0(species, '_gene_ensembl'))
  gos <- getBM(gene.list, attributes = c("ensembl_gene_id", 
                                         "external_gene_name"), filters = c("ensembl_gene_id"), 
               mart = mart.obj)
  gl <- gos[match(gene.list, gos[, 1]), 2]
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original) {
    gl[is.na(gl)] <- gene.list[is.na(gl)]
  }
  return(gl)
}



  # Load DE genes from Seurat -----------------------------------------------

jsuffix <- "seurat"

de.inf <- "/home/jyeung/data/from_rstudioserver/zebrafish/diff_exprs_Chloe_seurat.full.rds"
de.exprs <- readRDS(de.inf)

pmin <- 0.01
fcmin <- ""

# # load pseudobulk 
# jsuffix <- "cpmnorm"
# inf.bulk <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM_", jsuffix, ".2019-12-10.rds")
# # inf.bulk <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"
# wkm.pbulk <- readRDS(inf.bulk)



# Load TSS signal?  -------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

jmark <- jmarks[[1]]

# bad.ctype <- c("")
# bad.ctype.grep <- "eryth|Unknown"
bad.ctype.grep <- ""
btypes.keep <- c("tss")
# btypes.keep <- c("enh")

bad.ctype.grep.str <- paste(strsplit(bad.ctype.grep, split = "\\|")[[1]], collapse = "_")
btypes.keep.str <- paste(btypes.keep, collapse = "_")
outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis.", jsuffix)
dir.create(outdir)
for (jmark in jmarks){
  if (bad.ctype.grep == ""){
    pdfout <- file.path(outdir, paste0("pseudobulks.", btypes.keep.str, ".AllCtypes.", jmark, ".pdf"))
  }  else {
    pdfout <- file.path(outdir, paste0("pseudobulks.", btypes.keep.str, ".filter_", bad.ctype.grep.str, ".", jmark, ".pdf"))
  }
  
  if (file.exists(pdfout)){
    print(paste("pdfout exists, skipping:", pdfout))
    next
  }
  
  pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)
  
  
  jprefix <- "/home/jyeung/hub_oudenaarden"
  
  inf.tss <- file.path(jprefix, paste0("jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.winsize_10000/PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
  inf.enh <- file.path(jprefix, paste0("jyeung/data/zebrafish_scchic/count_tables_all/count_tables.ConservedEnhancers.winsize_10000/PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
  inf.genebody <- file.path(jprefix, paste0("jyeung/data/zebrafish_scchic/count_tables_all/count_tables.genebody/PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.genebody.csv"))
  
  inf.clstr <- file.path(jprefix, paste0("jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt"))
  
  mat.tss <- ReadMatTSSFormat(inf.tss, as.sparse=TRUE)
  mat.enh <- ReadMatTSSFormat(inf.enh, as.sparse=TRUE, add.coord = TRUE)
  rownames(mat.enh) <- paste("chr", rownames(mat.enh), sep = "")
  mat.gb.full <- ReadMatTSSFormat(inf.genebody, as.sparse=TRUE, add.coord = FALSE)
  
  mat.tss <- CollapseRowsByGene(mat.tss, as.long=FALSE, track.kept.gene = TRUE)
  mat.gb <- CollapseRowsByGene(mat.gb.full, as.long=FALSE, track.kept.gene = TRUE)
  
  
  clstrs <- fread(inf.clstr)
  
  
  # Assign enhancers to gene ------------------------------------------------
  
  
  # load bins
  # terms.mat.tmp <- tm.result$terms
  # rownames(terms.mat.tmp) <- rownames(terms.mat.tmp)
  
  jwin <- 50000L
  jchromos <- paste("chr", seq(25), sep = "")
  inf.annot <- file.path(jprefix, paste0("jyeung/data/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed"))
  assertthat::assert_that(file.exists(inf.annot))
  
  jcoords.vec.full <- sapply(rownames(mat.enh), function(x) strsplit(x, ";")[[1]][[1]], simplify = TRUE)
  
  # keep only uniques?
  mat.enh.filt <- mat.enh[!duplicated(jcoords.vec.full), ]
  rownames(mat.enh.filt) <- sapply(rownames(mat.enh.filt), function(x) strsplit(x, ";")[[1]][[1]])
  
  annot.out <- AnnotateCoordsFromList(coords.vec = rownames(mat.enh.filt), inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
  
  # assign region to gene
  jhash <- hash(annot.out$out2.df.closest$region_coord, annot.out$out2.df.closest$SYMBOL)
  
  mat.enh.filt[1:5, 1:5]
  
  rnames.new <- sapply(rownames(mat.enh.filt), function(x) paste(x, jhash[[x]], sep = ";"))
  rownames(mat.enh.filt) <- rnames.new
  
  # remove rows that end with ; (no gene names)
  mat.enh.filt <- mat.enh.filt[!grepl(";$", rownames(mat.enh.filt)), ]
  
  mat.lst <- list(tss = mat.tss, enh = mat.enh.filt, gb = mat.gb)
  
  btypes.remove <- names(mat.lst)[!names(mat.lst) %in% btypes.keep]
  
  print(paste("Removing btypes:"))
  print(btypes.remove)
  
  # remove btypes
  for (btype in btypes.remove){
    mat.lst[[btype]] <- NULL
  }
  
  # Create pseudobulks  -----------------------------------------------------
  
  # bad.ctype <- c("Unknown", "eryth")
  # bad.ctype <- c("Unknown")
  
  # # keep ;1 ? 
  # rnames.keep <- grepl(";1$", rownames(mat))
  # mat.filt <- mat[rnames.keep, ]
  
  cnames.keep.lst <- lapply(split(x = clstrs, f = clstrs$cluster), function(jdat) jdat$cell)
  
  mat.pbulk.lst <- lapply(mat.lst, function(x){
    mat.pbulk <- SumAcrossClusters(x, cnames.keep.lst) %>%
      bind_rows() %>%
      as.data.frame()
    rownames(mat.pbulk) <- rownames(x)
    
    meta <- data.table(ctype = colnames(mat.pbulk), stringsAsFactors = FALSE)
    dds <- DESeqDataSetFromMatrix(mat.pbulk, meta, design = ~ 1)
    dds.vst <- as.matrix(assay(vst(dds, fitType = "local")))
    return(dds.vst)
  })
  
  mat.pbulk.lst <- lapply(mat.pbulk.lst, function(x){
    # cnames.keep.bool <- !colnames(x) %in% bad.ctype
    if (bad.ctype.grep == ""){
      cnames.keep.bool <- rep(TRUE, ncol(x))
    } else {
      cnames.keep.bool <- !grepl(bad.ctype.grep, colnames(x))
    }
    return(x[, cnames.keep.bool])
  })
  
  dat.sum.long <- lapply(names(mat.pbulk.lst), function(btype){
    print(btype)
    jmat <- mat.pbulk.lst[[btype]]
    jdat.long <- melt(data.frame(genefull = rownames(jmat), jmat, stringsAsFactors = FALSE), variable.name = "ctype", value.name = "logexprs") %>%
      rowwise() %>%
      mutate(gene = strsplit(as.character(genefull), ";")[[1]][[2]]) %>%
      group_by(gene) %>%
      mutate(logFC = logexprs - mean(logexprs),
             zscore = logFC / sd(logexprs),
             btype = btype)
    return(jdat.long)
  }) %>%
    bind_rows()
  
  
  
  # Get ctype specific genes ------------------------------------------------
  
  # pick randomly 100 genes
  
  # get eryth-specific genes
  pvalmin <- 1e-3
  logfcmin <- 1
  
  print(unique(de.exprs$cluster))
  
  jctype <- "monocytes"
  jsub <- subset(de.exprs, cluster == jctype & p_val_adj <= pvalmin & avg_logFC > logfcmin)
  
  print(dim(jsub))
  
  gkeep <- sapply(jsub$gene, function(x) strsplit(x, "-")[[1]][[2]])
  
  gkeep <- sample(x = unique(dat.sum.long$gene), size = 100)
  
  m <- ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = ctype, y = zscore)) + geom_boxplot() + 
    geom_jitter(width = 0.15) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    facet_wrap(~btype) + ggtitle("Random genes:", jmark)
  
  
  # Define gene sets --------------------------------------------------------
  
  gbefore <- unique(as.character(wkm.pbulk$gene))
  gafter <- EnsemblGene2Gene.ZF(gene.list = gbefore, return.original = TRUE)
  ghash <- hash(gbefore, gafter)
  
  jctypes <- unique(wkm.pbulk$celltype)
  names(jctypes) <- jctypes
  
  gsets <- lapply(jctypes, function(jctype){
    jsub <- subset(wkm.pbulk, celltype == jctype & zscore > 1.5)
    jgenes <- as.character(jsub$gene)
    jgenes.gnames <- sapply(jgenes, function(x) ghash[[x]], USE.NAMES = FALSE)
    return(jgenes.gnames)
  })
  
  
  # Check zscore distribution for each gene set and each celltype -----------
  
  gset.names <- names(gsets)
  names(gset.names) <- gset.names
  
  gset.all <- unique(unlist(gsets))
  
  dat.sum.list.gsets <- lapply(gset.names, function(gset.name){
    print(gset.name)
    gset <- gsets[[gset.name]]
    dat.sum.long <- dat.sum.long %>%
      rowwise() %>%
      mutate(geneset = ifelse(gene %in% gset, gset.name, NA))
    return(subset(dat.sum.long, !is.na(geneset)))
  })
  
  # add the NAgroup group
  dat.sum.long.bg <- dat.sum.long %>%
    rowwise() %>%
    mutate(geneset = ifelse(gene %in% gset.all, "zBackground", NA)) %>%
    filter(!is.na(geneset))
  
  dat.sum.list.gsets$zBackground <- dat.sum.long.bg
  
  dat.sum.long.gsets <- bind_rows(dat.sum.list.gsets)
  
  
  m <- ggplot(dat.sum.long.gsets, aes(x = logexprs, group = btype, fill = btype)) + geom_density(alpha = 0.25) + 
    facet_grid(ctype ~ geneset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + 
    geom_vline(xintercept = 0)
  print(m)
  m <- ggplot(dat.sum.long.gsets, aes(x = logFC, group = btype, fill = btype)) + geom_density(alpha = 0.25) + 
    facet_grid(ctype ~ geneset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + 
    geom_vline(xintercept = 0)
  print(m)
  
  m <- ggplot(dat.sum.long.gsets, aes(x = zscore, group = btype, fill = btype)) + geom_density(alpha = 0.25) + 
    facet_grid(ctype ~ geneset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark) + 
    geom_vline(xintercept = 0)
  print(m)
  
  
  # random genes?
  gkeep <- sample(x = unique(dat.sum.long$gene), size = 100, replace = FALSE)
  
  m <- ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = zscore, fill = btype)) + geom_density(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark, "random genes") + 
    geom_vline(xintercept = 0)
  print(m)
  
  m <- ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = logFC, fill = btype)) + geom_density(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
    ggtitle(jmark, "random genes") + 
    geom_vline(xintercept = 0)
  print(m)
  
  m <- ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = logexprs, fill = btype)) + geom_density(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())   + 
    ggtitle(jmark, "random genes") + 
    geom_vline(xintercept = 0)
  print(m)
  
  m <- ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = ctype, y = zscore)) + geom_boxplot() + 
    geom_jitter(width = 0.15) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    ggtitle(jmark, "random genes") + 
    geom_vline(xintercept = 0)
  print(m)
  
  dev.off()
  
}




# ggplot(dat.sum.long.gsets %>% filter(ctype == "eryth" & geneset == "erythrocytes"), aes(x = zscore, fill = btype)) + 
#   geom_density(alpha = 0.25) + 
#   facet_grid(ctype ~ geneset) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 


# 
# dat.sum.long.lst <- lapply(split(dat.sum.long, f = dat.sum.long$ctype), function(jdat){
#   jdat$gene
# })
# 
# 
# # Try from pseudbulk ------------------------------------------------------
# 
# print(unique(wkm.pbulk$celltype))
# 
# ggplot(wkm.pbulk, aes(x = exprs, fill = celltype, group = celltype)) + geom_density() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~celltype) + geom_vline(xintercept = 5)
# 
# ggplot(wkm.pbulk, aes(x = zscore, fill = celltype, group = celltype)) + geom_density() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~celltype)
# 
# jctype <- "HSPCs"
# jsub <- wkm.pbulk %>%
#   group_by(gene) %>%
#   filter(exprs[which(celltype == jctype)] > 5 & mean(exprs[which(exprs != jctype)]) < 5)
# 
# jsub <- subset(wkm.pbulk, celltype == jctype & zscore > 1)
# 
# gkeep <- unique(jsub$gene)
# 
# ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = ctype, y = zscore)) + geom_boxplot() + 
#   geom_jitter(width = 0.15) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


# # find lymph-specific regions??
# 
# gkeep <- (subset(dat.sum.long, ctype == "lymph") %>%
#   ungroup() %>%
#   arrange(desc(zscore)))$gene[1:100]
# 
# ggplot(dat.sum.long %>% filter(gene %in% gkeep), aes(x = ctype, y = zscore)) + geom_boxplot() + 
#   geom_jitter(width = 0.15) + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# # ggplot(wkm.pbulk %>% filter(gene %in% gkeep), aes(x = zscore, fill = celltype, group = celltype)) + geom_density() + 
# #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
# #   facet_wrap(~celltype)
# 
# 
# 

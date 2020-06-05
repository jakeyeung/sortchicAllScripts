# Jake Yeung
# Date of Creation: 2020-06-01
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/7-calculate_high_dimensional_fold_changes.R
# Summarize fold changes: find the different movements of genes in a genome-wide manner? 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(DropletUtils)

library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

hubprefix <- "/home/jyeung/hub_oudenaarden"
jorg <- "org.Mm.eg.db"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/first_transcript_tss/gene_tss_winsize.10000.first_transcript.bed")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# Constants ---------------------------------------------------------------

    
# indir.bins <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.count_tables_bl_filt")
# indir.bins <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.first_transcript")
indir.bins.old <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow")
indir.bins <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2")
inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.stringentDE.forAvO/filtered_cells_for_pseudobulk.2020-05-31.rds")

assertthat::assert_that(dir.exists(indir.bins))


# Load meta ---------------------------------------------------------------


dat.meta <- readRDS(inf.meta)


# Load bins matrix  -------------------------------------------------------

mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  topics.keep <- dat.meta[[jmark]]
  # infs.csv <- list.files(path = indir.bins, pattern = paste0(jmark, ".*.csv"), full.names = TRUE)
  inf <- paste0(jmark, ".mapq_40.SlidingWindow_10000.blfiltered.csv")
  inf.mat <- file.path(indir.bins, inf)
  # ReadMatTSSFormat(inf.mat)
  mat <- ReadMatSlideWinFormat(inf.mat, add.chromo = FALSE)
  return(mat)
})

# load pseudobulks 
mats.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  topics.keep <- dat.meta[[jmark]]
  # infs.csv <- list.files(path = indir.bins, pattern = paste0(jmark, ".*.csv"), full.names = TRUE)
  inf <- paste0(jmark, ".mapq_40.SlidingWindow_10000.blfiltered.csv")
  inf.mat <- file.path(indir.bins, inf)
  # ReadMatTSSFormat(inf.mat)
  mat <- ReadMatSlideWinFormat(inf.mat, add.chromo = FALSE)
  return(mat)
})


# Check periodicities in the odds and evens  ------------------------------


jout <- table(mats.lst$H3K4me1[, 1])

plot(log2(jout))

# Select TSS to annotate bins later ---------------------------------------

# we did 10kb for zf, use H3K4me3 as reference to choose which TSS for gene
jwinsize <- 10000L
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS"

ref.mark <- "H3K4me3"
fname <- paste0(ref.mark, ".countTableTSS.mapq_40.TSS_", jwinsize, ".blfiltered.csv")
inf <- file.path(indir, fname)
tss.mat.ref <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = TRUE)
tss.mat.ref.filt <- CollapseRowsByGene(tss.mat.ref, as.long = FALSE, track.kept.gene = TRUE)
tss.keep.fromref <- rownames(tss.mat.ref.filt)

counts.lst <- apply(tss.mat.ref.filt, 2, table)



tss.coords <- sapply(tss.keep.fromref, function(x) strsplit(x, ";")[[1]][[1]])
# tss.names <- sapply(tss.keep.fromref, function(x) strsplit(x, ";")[[1]][[2]])

tss.dat <- data.frame(seqnames = sapply(tss.coords, GetChromo), 
                      start = sapply(tss.coords, GetStart), 
                      end = sapply(tss.coords, GetEnd), 
                      tssname = tss.keep.fromref,
                      stringsAsFactors = FALSE)


# Get pseudobulks  --------------------------------------------------------

# create cells into topics

cells.clusters.lst <- lapply(dat.meta, function(jdat){
  jsplit <- lapply(split(x = jdat, f = jdat$cluster.new), function(x){
    return(x$cell)
  })
})

pmat.lst <- lapply(jmarks, function(jmark){
  jmat <- mats.lst[[jmark]]
  cnames.keep.lst <- cells.clusters.lst[[jmark]]
  pmat <- SumAcrossClusters(jmat, cnames.keep.lst) %>%
    as.data.frame()
})

print(lapply(pmat.lst, dim))

# take common bins
bins.common <- Reduce(f = intersect, x = lapply(pmat.lst, rownames))

pmat.filt.lst <- lapply(pmat.lst, function(pmat){
  pmat[bins.common, ]
})

# Annotate bins  ----------------------------------------------------------

annot.bins <- AnnotateCoordsFromList.GeneWise(coords.vec = bins.common, inf.tss = tss.dat, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = jorg, chromos.keep = jchromos)

# Downsample  -------------------------------------------------------------

pmat.ds.lst <- lapply(pmat.filt.lst, function(pmat){
  # downsample 
  jprop <- min(colSums(pmat)) / colSums(pmat)
  rnames <- rownames(pmat)
  pmat.ds <- DropletUtils::downsampleMatrix(as.matrix(pmat), prop = jprop)
  rownames(pmat.ds) <- rnames
  return(pmat.ds)
})

# check downsampling
lapply(pmat.ds.lst, colSums)

# annotate bins

head(annot.bins$out2.df)

g2b <- subset(annot.bins$out2.df, select = c(region_coord, gene))
genes.uniq <- g2b$gene
g2e <- JFuncs::Gene2Ensembl.ZF(gene.list = genes.uniq, species = "mmusculus")
g2e.hash <- hash(genes.uniq, g2e)
b2g.hash <- hash(g2b$region_coord, g2b$gene)
b2e.hash <- hash(g2b$region_coord, sapply(g2b$gene, function(g) AssignHash(g, g2e.hash, null.fill = g)))

pmat.ds.annot.lst <- lapply(pmat.ds.lst, function(pmat){
  pmat.dat <- data.frame(bin = rownames(pmat), 
                         gene = sapply(rownames(pmat), function(bin) AssignHash(bin, b2g.hash, null.fill = NA)), 
                         ens = sapply(rownames(pmat), function(bin) AssignHash(bin, b2e.hash, null.fill = NA)),
                         pmat, 
                         stringsAsFactors = FALSE)
})


# Do differential analysis  -----------------------------------------------

# these need to be done in log and not linear because linear has a mean-variance relationship. 
# Also then displacement is directly interpretable as fold change
# pmat.ds.long.lst <- lapply(pmat.ds.annot.lst, function(pmat.ds){
pmat.ds.long.lst <- lapply(jmarks, function(jmark){
  pmat.ds <- pmat.ds.annot.lst[[jmark]]
  jlong <- pmat.ds %>%
    data.table::melt(data = ., id.vars = c("bin", "gene", "ens"), variable.name = "pseudobulk", value.name = "counts") %>%
    group_by(bin) %>%
    mutate(log2counts = log2(counts  + 1),
           log2fc = log2counts - mean(log2counts),
           fczscore = log2fc / sd(log2counts))
  jlong$mark <- jmark
  return(jlong)
}) %>%
  bind_rows()

# check counts
# plot(density(subset(pmat.ds.long.lst, mark == "H3K4me3")$counts))
# plot(density(log2(subset(pmat.ds.long.lst, mark == "H3K4me3")$counts + 1)))

for (jmark in jmarks){
  plot(density(subset(pmat.ds.long.lst, mark == jmark)$log2counts), main = jmark)
}

# make logcounts matrix and plot

ctypes <- as.character(unique(pmat.ds.long.lst$pseudobulk))
names(ctypes) <- ctypes
pmat.wide.byctype <- lapply(ctypes, function(ctype){
  jsub <- subset(pmat.ds.long.lst, pseudobulk == ctype)
  jsub.wide <- reshape2::dcast(data = jsub, formula = "bin ~ mark", value.var = "counts")
  jsub.wide$pseudobulk <- ctype
  return(jsub.wide)
})


# Plot 2D plots -----------------------------------------------------------

mlst <- lapply(ctypes, function(ctype){
  m <- ggplot(pmat.wide.byctype[[ctype]], aes(x = H3K4me3, y = H3K27me3)) + 
    geom_point(alpha = 0.2) + 
    geom_density_2d() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggtitle(ctype)
  return(m)
})

print(mlst)

# Load gene sets ----------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)








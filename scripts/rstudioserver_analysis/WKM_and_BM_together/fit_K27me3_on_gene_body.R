# Jake Yeung
# Date of Creation: 2020-08-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/fit_K27me3_on_gene_body.R
# Fit on gene body

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Load mats  --------------------------------------------------------------

jmark.test <- "H3K27me3"
hubprefix <- "/home/jyeung/hub_oudenaarden"
inf <- paste0(hubprefix, "/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40/countTablesAndRZr1only.NewFilters/", jmark.test, "-BM_AllMerged.merged_by_clusters_with_NAs.count_table.TSS_TES.csv")

mat <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = TRUE)



# Load annots -------------------------------------------------------------

RenameClusterBM <- function(clstr.orig, bm.rename){
  # clstr.orig <- "Bcells-Cd83_topic10"
  clstr.new <- paste0("z", clstr.orig)
  for (cname in names(bm.rename)){
    if (startsWith(clstr.orig, prefix = cname)){
      clstr.new <- bm.rename[[cname]]
    } else{
    }
  }
  return(clstr.new)
}

bm.rename <- as.list(hash(c("Bcells", "Eryth", "HSCs", "Neutrophils"), c("lymph", "eryth", "HSPCs", "granu")))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
dat.annot.lst.BM <- lapply(jmarks, function(jmark){
  inf.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
  assertthat::assert_that(file.exists(inf.annot))
  load(inf.annot, v=T)
  dat.umap.glm.fillNAs <- subset(dat.umap.glm.fillNAs, !is.na(cluster))
  dat.umap.glm.fillNAs$cluster <- sapply(dat.umap.glm.fillNAs$cluster, RenameClusterBM, bm.rename)
  return(dat.umap.glm.fillNAs)
})



# Get example genes  ------------------------------------------------------

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)

head(de.ens.sorted.stringent$Bcell)
e2g.hash <- hash::invert(g2e.hash)
evec <- as.character(de.ens.sorted.stringent$Erythroblast)
gvec <- sapply(evec, function(x) AssignHash(x, e2g.hash, null.fill = NA))

genes.annot <- data.frame(ens = names(e2g.hash), gene = hash::values(e2g.hash), stringsAsFactors = FALSE)

g2r.annot <- data.frame(row = rownames(mat), gene = sapply(rownames(mat), function(x) strsplit(strsplit(x, "\\.\\.")[[1]][[2]], split = "_")[[1]][[1]]), stringsAsFactors = FALSE)

gre.annot <- left_join(genes.annot, g2r.annot)

# Plot a neutro gene  -----------------------------------------------------

# jgene <- "S100a9"
# jgene <- "Ebf1"
# jgene <- ""
# (jgene <- sample(x = gvec[!is.na(gvec)], size = 1))
# jgene <- "Hlf"

(jrow <- grep(jgene, rownames(mat), value = TRUE))

(jrow <- sample(rownames(mat), 1))
(jgene <- strsplit(jrow, "\\.\\.")[[1]][[2]])

jvec <- mat[grepl(jgene, rownames(mat)), ]

cellsums <- data.frame(cell = colnames(mat), ncuts = colSums(mat), stringsAsFactors = FALSE)

# plot fits

input.dat <- data.frame(genecounts = mat[jrow, ], cell = names(jvec)) %>%
  left_join(., cellsums) %>%
  left_join(., dat.annot.lst.BM[[jmark.test]])

ggplot(input.dat, aes(x = cluster, y = log2(genecounts / ncuts))) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle(jgene, jrow)


cnames.keep.lst <- split(dat.annot.lst.BM[[jmark.test]], 
                         f = (dat.annot.lst.BM[[jmark.test]]$cluster))

for (jname in names(cnames.keep.lst)){
  if (startsWith("^z", jname)){
    cnames.keep.lst[[jname]] <- NULL
  } else {
    cnames.keep.lst[[jname]] <- cnames.keep.lst[[jname]]$cell
  }
}

mat.sum <- SumAcrossClusters(mat, cnames.keep.lst)
mat.sum <- do.call(cbind, mat.sum)


mat.sum <- log2(sweep(mat.sum, MARGIN = 2, STATS = colSums(mat.sum), FUN = "/") * 10^6 + 1) 

mat.sum.long <- mat.sum %>%
  melt()
colnames(mat.sum.long) <- c("gene", "cluster", "logexprs")

mat.sum.long <- mat.sum.long %>%
  group_by(cluster) %>%
  mutate(zscore = scale(logexprs, center = TRUE, scale = TRUE))

ggplot(mat.sum.long, aes(x = cluster, y = zscore)) +
  geom_point()  + 
  geom_boxplot() + 
  theme_bw()  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

 
 
jgene <- "Tal1"
(jrow <- grep(jgene, rownames(mat), value = TRUE))

# jrow <- sample(rownames(mat), size = 1)
# (jgene <- strsplit(jrow, "\\.\\.")[[1]][[2]])

print(names(de.ens.sorted.stringent))

jgset <- "Bcell"
jgset <- "Neutrophil"
jgset <- "Erythroblast"

evec <- as.character(de.ens.sorted.stringent[[jgset]])
# gvec <- sapply(evec, function(x) AssignHash(x, e2g.hash, null.fill = NA))
rowsvec <- subset(gre.annot %>% filter(!is.na(row)), ens %in% evec)$row

ggplot(mat.sum.long %>% filter(gene %in% rowsvec), aes(x = cluster, y = zscore)) +
  geom_point()  + 
  geom_boxplot() + 
  theme_bw()  + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  ggtitle(jgset)
  


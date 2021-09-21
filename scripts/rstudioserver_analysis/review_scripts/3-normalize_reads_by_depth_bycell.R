# Jake Yeung
# Date of Creation: 2021-06-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/3-normalize_reads_by_depth_bycell.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(hash)

hubprefix <- "/home/jyeung/hub_oudenaarden"
cnames <- c("fname1", "fname2", "nlines1", "nlines2")



# Get zeller meta  --------------------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K562_clean_by_G1filt.fracnonzerofilt/K562_cleaned_by_G1filt.H3K27me3.nonzerofracthres_2.bsize_50000.txt")
dat.meta <- fread(inf.meta)
good.cells <- dat.meta$cell



# Get good table ----------------------------------------------------------

inrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-09.rds")
dat.all <- readRDS(inrds, file = inrds)



# Normalize kaya-okur  ----------------------------------------------------

indir.ko <- file.path(hubprefix, "jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/count_before_after_filtering_hiddendomains_filt")
infs.ko <- list.files(indir.ko, pattern = "*.txt", full.names = TRUE)

dat.ko.dup <- lapply(infs.ko, function(inf){
  fread(inf, header = FALSE, col.names = cnames)
}) %>%
  bind_rows()

# get norm factor 
# avg reads per cell dat.ko.dup

# ntotal by cells 
dat.ko.jfac <- dat.ko.dup %>%
  group_by(fname1) %>%
  summarise(ntotal = sum(nlines1)) %>%
  mutate(jset = "KayaOkur") %>%
  dplyr::rename(cell = fname1) %>%
  mutate(cell = basename(cell),
         cell = gsub("bed.gz", "dedup.bed", cell))
  

# Normalize Wu et al  -----------------------------------------------------

# hESC
inf.hesc <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/GSM4780540_K27me3_H1_r1.fragments.HG38.tsv")
dat.hesc <- fread(inf.hesc, col.names = c("chromo", "start", "end", "cell", "dupcount")) %>%
  ungroup() %>%
  mutate(count = 1)

dat.hesc.sum <- dat.hesc %>%
  group_by(cell) %>%
  summarise(dupcount = sum(dupcount),
            count = sum(count))

ggplot(dat.hesc.sum, aes(x = count / dupcount)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.hesc.jfac <- dat.hesc %>%
  group_by(cell) %>%
  summarise(ntotal = sum(dupcount)) %>%
  mutate(jset = "Wu2021_hESC") 


# Normalize Bartosovic ----------------------------------------------------

inf.bart <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/GSE157637/GSE157637_H3K27me3_cell_lines_seurat_object.Rds")
dat.bart.seurat <- readRDS(inf.bart)
cell.ids <- sapply(as.character(Idents(dat.bart.seurat)), function(x) ClipLast(x, jsep = "_"))
cell.ids.hash <- hash::hash(names(cell.ids), as.character(cell.ids))

cell.ids.hash.names <- hash::hash(names(Idents(dat.bart.seurat)), as.character(cell.ids))

ctypes <- c("mESC", "3T3", "Oli-neu")
ctype2set <- hash::hash(ctypes, interaction("Bart2021", ctypes, sep = "_"))

inf.bart.frags <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/GSE157637/GSE157637_fragments.tsv")
dat.bart.frags <- fread(inf.bart.frags, header = FALSE, col.names = c("chromo", "start", "end", "cell", "dupcounts"))

dat.bart.jfac <- dat.bart.frags %>%
  group_by(cell) %>%
  summarise(ntotal = sum(dupcounts)) %>%
  rowwise() %>%
  mutate(ctype = AssignHash(x = cell, jhash = cell.ids.hash.names, null.fill = cell)) %>%
  filter(ctype %in% unique(cell.ids)) %>%
  mutate(jset = AssignHash(x = ctype, jhash = ctype2set, null.fill = cell)) %>%
  dplyr::select(-ctype)
  

# Normalize zeller  -------------------------------------------------------

inf.zeller.nodedup <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters_NoDedup.check_dedup_vs_nodedup/K562_AllMerged_H3K27me3.merged.sorted.tagged.countTable.binsize_50000.nodedup.csv")
inf.zeller.dedup <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters_NoDedup.check_dedup_vs_nodedup/K562_AllMerged_H3K27me3.merged.sorted.tagged.countTable.binsize_50000.dedup.csv")
mat.zeller.nodedup <- ReadMatSlideWinFormat(inf.zeller.nodedup)
mat.zeller.dedup <- ReadMatSlideWinFormat(inf.zeller.dedup)

dat.zeller.jfac <- data.frame(cell = colnames(mat.zeller.nodedup), ntotal = colSums(mat.zeller.nodedup), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(jset = "Zeller")

dat.zeller.jfac.dedup <- data.frame(cell = colnames(mat.zeller.dedup), ntotal = colSums(mat.zeller.dedup), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(jset = "Zeller.dedup")

dat.zeller.jfac.merge <- left_join(dat.zeller.jfac, dat.zeller.jfac.dedup, by = "cell")

ggplot(dat.zeller.jfac.merge, aes(x = ntotal.y / ntotal.x)) + geom_density()  + 
  xlab("Number of Unique Reads / Number of Total Reads")




# Make hash  --------------------------------------------------------------

dat.all.jfac <- rbind(dat.ko.jfac, 
                      dat.hesc.jfac, 
                      dat.bart.jfac,
                      dat.zeller.jfac)


# Fix dat.all  ------------------------------------------------------------

cells.keep <- subset(dat.all, jset == "Zeller")$cell
dat.all.tmp <- subset(dat.all, jset != "Zeller")
dat.all.toadd <- dat.zeller.jfac.dedup %>%
  dplyr::rename(nreads = ntotal) %>%
  mutate(jset = "Zeller") %>% 
  filter(cell %in% cells.keep)
dat.all.fixed <- rbind(dat.all.tmp, dat.all.toadd) 

outrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.", Sys.Date(), ".Zeller_dedup_fixed.rds")
saveRDS(dat.all.fixed, file = outrds)


# Replot  -----------------------------------------------------------------

m1 <- ggplot(dat.all.fixed, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = log2(nreads))) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("log2(unique reads)") + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# dat.all <- 

dat.all.merge <- left_join(dat.all.fixed, dat.all.jfac, by = c("cell", "jset")) %>%
  filter(jset != "Ku_WBCs")

ggplot(dat.all.merge, aes(forcats::fct_reorder(jset, nreads/ntotal, median, .desc = TRUE), y = log2(nreads / ntotal))) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("log2(unique reads / ntotal)") + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m.cellspecificnorm <- ggplot(dat.all.merge, aes(x = ntotal, y = nreads, color = jset)) + 
  geom_point() + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(slope = 1, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


outrds.merged <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.", Sys.Date(), ".Zeller_dedup_fixed.merged_with_cellspec_norm.rds")
saveRDS(dat.all.merge, file = outrds.merged)




# 
# # Replot noramlized  ------------------------------------------------------
# 
# # dat.all.norm <- left_join(dat.all, dat.all.jfac) %>%
# #   filter(!is.na(jfac))
# 
# # plot the correction factor 
# cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# m.factors <- ggplot(dat.all.jfac, aes(x = ntotal, y = ncells, color = jset)) + 
#   scale_color_manual(values = cbPalette) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_x_log10() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# m.factors2 <- ggplot(dat.all.jfac, aes(y = ntotal / ncells, x = jset)) + 
#   geom_col() + 
#   theme_bw() + 
#   ylab("total_reads per cell") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# m2 <- ggplot(dat.all.norm, aes(x = forcats::fct_reorder(jset, nreads / jfac, median, .desc = TRUE), y = log2(nreads / jfac))) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   ylab("log2(unique_reads / total_reads per cell") + 
#   xlab("") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# m3 <- ggplot(dat.all.norm, aes(x = forcats::fct_reorder(jset, nreads / jfac, median, .desc = TRUE), y = log2(nreads / (jfac / 1000)))) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   ylab("log2(unique_reads / 1000_total_reads per cell") + 
#   xlab("") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# # m4 <- ggplot(dat.all.norm, aes(x = nreads, y = ntotal))) + 
# #   geom_boxplot() + 
# #   theme_bw() + 
# #   ylab("log2(unique_reads / 1000_total_reads per cell") + 
# #   xlab("") + 
# #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 




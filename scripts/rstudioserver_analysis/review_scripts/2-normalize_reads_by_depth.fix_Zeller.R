# Jake Yeung
# Date of Creation: 2021-06-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/2-normalize_reads_by_depth.R
# Try to normalize by sequencing depth 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"
cnames <- c("fname1", "fname2", "nlines1", "nlines2")



# Get zeller meta  --------------------------------------------------------

inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K562_clean_by_G1filt.fracnonzerofilt/K562_cleaned_by_G1filt.H3K27me3.nonzerofracthres_2.bsize_50000.txt")
dat.meta <- fread(inf.meta)
good.cells <- dat.meta$cell



# Get good table ----------------------------------------------------------

# inrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-09.rds")
inrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-12.Zeller_dedup_fixed.rds")
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
ntotal.ko <- sum(dat.ko.dup$nlines1)
ncells.ko <- length(subset(dat.all, jset == "KayaOkur")$cell)
jfac.ko <- ntotal.ko / ncells.ko

dat.ko.jfac <- data.frame(ntotal = ntotal.ko, 
                          ncells = ncells.ko, 
                          jfac = jfac.ko, 
                          jset = "KayaOkur", 
                          stringsAsFactors = FALSE)

# Normalize Wu et al  -----------------------------------------------------

# hESC
inf.hesc <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/GSM4780540_K27me3_H1_r1.fragments.HG38.tsv")
dat.hesc <- fread(inf.hesc, col.names = c("chromo", "start", "end", "cell", "dupcount")) %>%
  ungroup() %>%
  mutate(count = 1)

dat.hesc.jfac <- dat.hesc %>%
  ungroup() %>%
  summarise(ntotal = sum(dupcount), 
            ncells = length(subset(dat.all, jset == "Wu2021_hESC")$cell)) %>%
  mutate(jfac = ntotal / ncells, 
         jset = "Wu2021_hESC")


# Normalize Bartosovic ----------------------------------------------------


inf.bart <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/GSE157637/GSE157637_H3K27me3_cell_lines_seurat_object.Rds")
dat.bart.seurat <- readRDS(inf.bart)
cell.ids <- sapply(as.character(Idents(dat.bart.seurat)), function(x) ClipLast(x, jsep = "_"))
cell.ids.hash <- hash::hash(names(cell.ids), as.character(cell.ids))

inf.bart.frags <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/GSE157637/GSE157637_fragments.tsv")
dat.bart.frags <- fread(inf.bart.frags, header = FALSE, col.names = c("chromo", "start", "end", "cell", "dupcounts"))

nreads.total <- sum(dat.bart.frags$dupcounts)

# divide it up by number lf cells 
ncells <- as.vector(unlist(table(cell.ids)))

jfac.bart <- nreads.total / sum(ncells)

nreads.vec <- nreads.total * (ncells / sum(ncells))

dat.bart.jfac <- data.frame(ntotal = nreads.vec, 
                            ncells = ncells,
                            jfac = jfac.bart, 
                            jset = c("Bart2021_Oli-neu", "Bart2021_mESC", "Bart2021_3T3"), 
                            stringsAsFactors = FALSE)


# Normalize zeller  -------------------------------------------------------


# inf.zeller <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters_NoDedup/K562_AllMerged_H3K27me3.merged.sorted.tagged.countTable.binsize_50000.csv")
inf.zeller <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters_NoDedup.check_dedup_vs_nodedup/K562_AllMerged_H3K27me3.merged.sorted.tagged.countTable.binsize_50000.nodedup.csv")
# inf.zeller <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/countTablesAndRZr1only_ByBins.NewFilters_NoDedup/K562_AllMerged_H3K27me3.merged.sorted.tagged.countTable.binsize_50000.csv")
mat.zeller <- ReadMatSlideWinFormat(inf.zeller)

dat.zeller.jfac <- data.frame(ntotal = sum(mat.zeller), ncells = length(good.cells), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(jfac = ntotal / ncells, 
         jset = "Zeller")


# Make hash  --------------------------------------------------------------

library(hash)
dat.all.jfac <- rbind(dat.ko.jfac, 
                      dat.hesc.jfac, 
                      dat.bart.jfac,
                      dat.zeller.jfac)

# Replot  -----------------------------------------------------------------

m1 <- ggplot(dat.all, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = log2(nreads))) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("log2(unique reads)") + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Replot noramlized  ------------------------------------------------------

dat.all.norm <- left_join(dat.all, dat.all.jfac) %>%
  filter(!is.na(jfac))

# plot the correction factor 
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.factors <- ggplot(dat.all.jfac, aes(x = ntotal, y = ncells, color = jset)) + 
  scale_color_manual(values = cbPalette) + 
  geom_point() + 
  theme_bw() + 
  scale_x_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m.factors2 <- ggplot(dat.all.jfac, aes(y = ntotal / ncells, x = jset)) + 
  geom_col() + 
  theme_bw() + 
  ylab("total_reads per cell") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

m2 <- ggplot(dat.all.norm, aes(x = forcats::fct_reorder(jset, nreads / jfac, median, .desc = TRUE), y = log2(nreads / jfac))) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab("log2(unique_reads / total_reads per cell") + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

m3 <- ggplot(dat.all.norm, aes(x = forcats::fct_reorder(jset, nreads / jfac, median, .desc = TRUE), y = log2(nreads / (jfac / 1000)))) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab("log2(unique_reads / 1000_total_reads per cell") + 
  xlab("") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# m4 <- ggplot(dat.all.norm, aes(x = nreads, y = ntotal))) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   ylab("log2(unique_reads / 1000_total_reads per cell") + 
#   xlab("") + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/CUTnTAG_vs_Zeller_H3K27me3_dedup.", Sys.Date(), ".zeller_fixed.pdf")
pdf(file = outpdf, useDingbats = FALSE)

print(m1 + ggtitle("H3K27me3 sortChIC vs CUTnTAG"))
print(m.factors + ggtitle("H3K27me3 sortChIC vs CUTnTAG"))
print(m.factors2 + ggtitle("H3K27me3 sortChIC vs CUTnTAG"))
print(m2 + ggtitle("H3K27me3 sortChIC vs CUTnTAG"))
print(m3 + ggtitle("H3K27me3 sortChIC vs CUTnTAG"))

dev.off()


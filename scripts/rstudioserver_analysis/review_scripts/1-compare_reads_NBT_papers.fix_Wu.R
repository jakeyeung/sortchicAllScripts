# Jake Yeung
# Date of Creation: 2021-06-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/1-compare_reads_NBT_papers.fix_Wu.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Bart --------------------------------------------------------------


inf.bart <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/GSE157637/GSE157637_H3K27me3_cell_lines_seurat_object.Rds")
dat.bart.seurat <- readRDS(inf.bart)

cell.ids <- Idents(dat.bart.seurat)
count.mat <- dat.bart.seurat@assays$bins_5000@counts

dat.cellsum <- data.frame(cell = names(cell.ids), ctype = cell.ids, nreads = colSums(count.mat), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(ctype.clipped = ClipLast(as.character(ctype), jsep = "_"),
         paper = "Bart2021",
         jset = interaction(paper, ctype.clipped, sep = "_"))

ggplot(dat.cellsum, aes(x = jset, y = log2(nreads))) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Wu et al ----------------------------------------------------------------

# hESC
inf.hesc <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/GSM4780540_K27me3_H1_r1.fragments.HG38.tsv")
dat.hesc <- fread(inf.hesc, col.names = c("chromo", "start", "end", "cell", "dupcount")) %>%
  ungroup() %>%
  mutate(count = dupcount)

dat.hesc.sum <- dat.hesc %>%
  group_by(cell) %>%
  summarise(nreads = sum(count)) %>%
  mutate(jset = "Wu2021_hESC")

ggplot(dat.hesc.sum, aes(y = log2(nreads))) + 
  geom_boxplot() 

# Kaya okur  --------------------------------------------------------------

inrds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/dat_hdfilt_merge.KayaOkur_vs_Zeller.H3K27me3.2021-06-07.rds"
dat.hdfilt.merge <- readRDS(inrds) %>%
  dplyr::rename(nreads = cuts_in_chromo, 
                cell = samp,
                jset = set)


# Ku et al ----------------------------------------------------------------

inf.ku <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Ku_et_al_2021/GSE139857_scK27_4lanes_rc_at_79100peaks_width10k.txt")
dat.ku <- fread(inf.ku, header = FALSE)

dat.ku.sum <- data.frame(cell = paste("cell", seq(ncol(dat.ku)), sep = ""), nreads = colSums(dat.ku), stringsAsFactors = FALSE) %>%
  mutate(jset = "Ku_WBCs")

# All together ------------------------------------------------------------

dat.all <- rbind(dat.cellsum %>% dplyr::select(cell, nreads, jset), 
                 dat.hesc.sum %>% dplyr::select(cell, nreads, jset),
                 dat.ku.sum %>% dplyr::select(cell, nreads, jset),
                 dat.hdfilt.merge %>% dplyr::select(cell, nreads, jset))

ggplot(dat.all, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = log2(nreads))) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("H3K27me3 comparison across datasets") + 
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.all, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = log10(nreads))) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("H3K27me3 comparison across datasets") + 
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.all, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.all, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_y_log10() + 
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# save output
outrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.", Sys.Date(), ".fix_Wu.rds")
saveRDS(dat.all, file = outrds)

# Jake Yeung
# Date of Creation: 2020-10-26
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_merged/add_hoescht_signal_to_glmpca.check_from_bed.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

inf.hoescht <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.hoesch_staining/hoesch_on_K562_plates.rds"
dat.hoescht.all <- readRDS(inf.hoescht)
dat.hoescht <- subset(dat.hoescht.all, grepl(jmark.hoescht, x = experi))


# Load spikein data -------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.spike <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData")
load(inf.spike, v=T)

# Load data ---------------------------------------------------------------


jmark <- "H3K4me1"
jmark.hoescht <- "4me1.*.4.*_index$"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks.hoescht.prefix <- c("4me1", "4me3", "27me3", "9me3")
hoescht.suffix <- c(".*.4.*index$")
jmarks.hoescht <- paste(jmarks.hoescht.prefix, hoescht.suffix, sep = "")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562_cellcycle.hoescht.allmerged_frombed"
dir.create(outdir)


# Load countmatrix raw ----------------------------------------------------





# Filter good cells -------------------------------------------------------




# Add hoescht staining ----------------------------------------------------




# Plot --------------------------------------------------------------------

# 
# 
# 
# 
# for (i in seq(length(jmarks))){
#   
#   jmark <- jmarks[[i]]
#   jmark.hoescht <- jmarks.hoescht[[i]]
#   
#   outpdf <- file.path(outdir, paste0(jmark, "_cell_cycle_GLMPCA_with_hoescht2.pdf"))
#   pdf(outpdf, useDingbats = FALSE)
#   
#   
#   inf.glmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/glmpcaPois_K562_spikein/K562_count_tables_50000.", jmark, jprefix, jsuffix))
#   assertthat::assert_that(file.exists(inf.glmpca))
#   
#   load(inf.glmpca, v=T)
#   
#   
#   # Load hoescht annots -----------------------------------------------------
#   
#   
#   
#   # Plot UMAP  --------------------------------------------------------------
#   
#   jsettings <- umap.defaults
#   jsettings$n_neighbors <- 30
#   jsettings$min_dist <- 0.1
#   jsettings$random_state <- 123
#   
#   dat.umap <- DoUmapAndLouvain(topics.mat = glmpcaout$factors, jsettings = jsettings)
#   
#   dat.umap$dim1 <- glmpcaout$factors$dim1
#   dat.umap$dim2 <- glmpcaout$factors$dim2 
#   dat.umap$dim3 <- glmpcaout$factors$dim3
#   dat.umap$dim4 <- glmpcaout$factors$dim4
#   
#   dat.spikeins.mat <- scchicFuncs::AddCellCycleLabel.bydat(dat.spikeins.mat)
#   
#   totalcounts <- data.frame(cell = colnames(Y), totalcounts = colSums(Y), stringsAsFactors = FALSE)
#   dat.spikeins.mat <- left_join(dat.spikeins.mat, totalcounts, by = c("samp" = "cell")) %>%
#     filter(!is.na(totalcounts))
#   
#   
#   dat.umap.annot <- left_join(dat.umap, dat.spikeins.mat, by = c("cell" = "samp"))
#   
#   
#   
#   # add hoescht -------------------------------------------------------------
#   
#   dat.umap.annot <- dat.umap.annot %>%
#     rowwise() %>%
#     mutate(experi2 = strsplit(cell, "_")[[1]][[1]])
#   experi.str <- unique(dat.umap.annot$experi2)
#   assertthat::assert_that(length(experi.str) == 1)
#   
#   
#   # rename experi to match glmpca
#   dat.hoescht <- dat.hoescht %>%
#     rowwise() %>%
#     mutate(experi2 = experi.str)
#   
#   jdat <- left_join(dat.umap.annot, dat.hoescht, by = c("rowcoord" = "row.indx", "colcoord" = "col.indx", "experi2" = "experi2"))
#   
#   cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
#   
#   m.cc <- ggplot(dat.umap.annot, aes(x = dim1, y = dim2, color = cellcycle.str)) + 
#     geom_point()  + 
#     scale_color_manual(values = cbPalette) + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     ggtitle(jmark, "labeled by cellcycle")
#   
#   m.hoescht <- ggplot(jdat, aes(x = dim1, y = dim2, color = hoesch)) + 
#     geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_color_viridis_c() + 
#     ggtitle(jmark, "labeled by hoesch staining")
#   
#   m.spikeins <- ggplot(jdat, aes(x = dim1, y = dim2, color = log2(totalcounts / spikeincounts))) + 
#     geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     scale_color_viridis_c() + 
#     ggtitle(jmark, "colored by logratio  totalcounts / spikeincounts")
#   
#   m.spikeins.vs.hoescht <- ggplot(jdat, aes(x = hoesch, y = log2(totalcounts / spikeincounts), color = cellcycle)) + 
#     geom_point() + 
#     theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#     ggtitle(jmark, "colored by logratio  totalcounts / spikeincounts")
#   
#   print(m.cc)
#   print(m.hoescht)
#   print(m.spikeins)
#   print(m.spikeins.vs.hoescht)
#   
#   # ggplot(jdat, aes(x = dim2, y = dim3, color = hoesch)) + 
#   #   geom_point() + 
#   #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   #   scale_color_viridis_c() + 
#   #   ggtitle(jmark, "labeled by hoesch staining")
#   dev.off()
#   
# }
# 
# 


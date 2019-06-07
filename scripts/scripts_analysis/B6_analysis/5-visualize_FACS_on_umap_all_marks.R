# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/5-visualize_FACS_on_umap_all_marks.R
# All marks


library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(JFuncs)

source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load umaps  -------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

inmain <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.dats <- lapply(jmarks, function(jmark) file.path(inmain, paste0("dat_umap_long_with_louvain.", jmark, ".RData")))

dat.umap.long <- lapply(inf.dats, LoadUmap) %>%
  bind_rows()


# Plot umaps for sanity ---------------------------------------------------

PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K4me1"), xvar = "umap1", yvar = "umap2")
PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K4me3"), xvar = "umap1", yvar = "umap2")
PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K27me3"), xvar = "umap1", yvar = "umap2")
PlotXYNoColor(dat.umap.long %>% filter(mark == "H3K9me3"), xvar = "umap1", yvar = "umap2")


# Check the FACS data -----------------------------------------------------

dat.facs.filt <- lapply(jmarks, LoadFACSGetLoadings) %>%
  bind_rows()

# bind to umap long
dat.merge <- left_join(dat.umap.long, dat.facs.filt %>% dplyr::select(cell, loadings, mark)) %>%
  filter(!is.na(loadings))


# Save object for Chloe ---------------------------------------------------

head(dat.facs.filt)
outf <- "~/data/scchic/robjs/B6_objs/umap_and_facs_data_merged.Rdata"
if (!file.exists(outf)){
  save(dat.facs.filt, dat.umap.long, file = outf)
}


# Make sure loadings are in right region ----------------------------------

pdfout <- "/Users/yeung/data/scchic/pdfs/B6_figures/facs_on_umap/facs_on_umap.pdf"
pdf(pdfout, useDingbats = FALSE)
mlst <- lapply(jmarks, function(jmark){
  m <- PlotXYWithColor(dat.merge %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "loadings", jtitle = jmark, jsize = 4) 
  # print(m)
})
multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 4)
dev.off()

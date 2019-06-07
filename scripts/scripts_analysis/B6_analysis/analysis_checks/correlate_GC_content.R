# Jake Yeung
# Date of Creation: 2019-05-28
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/analysis_checks/correlate_GC_content.R
# Find gene rich and gene poor regions



# BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(GGally)
library(colorRamps)


source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


indir <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6"
inf.traj <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/B6_objs/traj_objs_all_marks.Rdata"

assertthat::assert_that(dir.exists(indir))
assertthat::assert_that(file.exists(inf.traj))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1){
  outdir <- args[[1]]
  assertthat::assert_that(dir.exists(outdir))
} else {
  outdir <- "/tmp"
}
pdfout <- file.path(outdir, "Fig4_SharpeningOverTraj.pdf")


# Constants ---------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"

# Load LDA objects --------------------------------------------------------

# indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

# inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)


# Make long for chromosome 15 ---------------------------------------------

jpseudo <- 0
jfac <- 10^6
jstr <- "chr15:"

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>%
  bind_rows()


# Do one trajectory -------------------------------------------------------

chunksize <- 2e7  # 20 MB
trajnames <- c("lymphoid", "granu", "eryth", "mega")
# jtraj <- "eryth"


# Add GC content  ---------------------------------------------------------

library(genoset)

coords <- unique(mat.sub.merge$coord)
gr.dat <- data.frame(seqnames = sapply(coords, GetChromo), start = sapply(coords, GetStart), end = sapply(coords, GetEnd))
gr <- GenomicRanges::makeGRangesFromDataFrame(gr.dat)
gr.gc <- calcGC(object = gr, bsgenome = BSgenome.Mmusculus.UCSC.mm10)

gr.gc.dat <- data.frame(chromo = sapply(names(gr.gc), GetChromo), 
                        start = as.numeric(sapply(names(gr.gc), GetStart)),
                        end = as.numeric(sapply(names(gr.gc), GetEnd)),
                        gc = gr.gc) %>%
  mutate(midpt = start + (end - start) / 2)

ggplot(gr.gc.dat %>% filter(midpt > 40e6 & midpt < 60e6), aes(x = midpt / 10^6, y = gc)) + geom_line() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Compare with chromosome  ------------------------------------------------

# save object

save(gr.gc.dat, file = "/Users/yeung/data/scchic/robjs/gr_gc_dat.RData")

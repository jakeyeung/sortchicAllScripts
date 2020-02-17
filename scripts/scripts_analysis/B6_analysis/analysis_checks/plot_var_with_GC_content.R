# Jake Yeung
# Date of Creation: 2020-01-14
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/analysis_checks/plot_var_with_GC_content.R
# description

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(JFuncs)
library(scchicFuncs)

library(GGally)

# source("scripts/Rfunctions/VariabilityFunctions.R")
# source("scripts/Rfunctions/Aux.R")
# source("scripts/Rfunctions/AuxB6.R")
# source("scripts/Rfunctions/PlotFunctions.R")
# 

# Constants ---------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"

# Load LDA objects --------------------------------------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)




# Make long for chromosome 15 ---------------------------------------------

jpseudo <- 0
jfac <- 10^6
jstr <- "chr15:"

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>%
  bind_rows()

# Layer on top GC/AT content  ---------------------------------------------

gc.inf <- "/Users/yeung/data/scchic/robjs/gr_gc_dat.RData"
load(gc.inf, v=T)

# add GC to data
gr.gc.dat <- gr.gc.dat %>%
  mutate(pos = start) %>%
  mutate(analysis = "gc")


# pdf(file = "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/correlate_with_gc.pdf", useDingbats = FALSE)

m.gc <- ggplot(data = gr.gc.dat %>% mutate(gc = ifelse(gc < quantile(gc, 0.01), quantile(gc, 0.01), gc)), aes(x = pos / 10^6, y = analysis, fill = gc)) + geom_tile() + theme_bw() + 
    theme_bw() + theme(aspect.ratio=0.05, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") +
    scale_fill_gradientn(colours = colorRamps::matlab.like(5), name = "GC") +
    ylab("") +
    xlab("Position (MB)")
print(m.gc)
print(m.slopes)

multiplot(m.gc, m.slopes, cols = 1)


# get variance for K4me3

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.impute.log <- log2(t(tm.result.lst$H3K4me3$topics %*% tm.result.lst$H3K4me3$terms))
dat.var <- CalculateVarAll(dat.impute.log = dat.impute.log, jchromos = jchromos)

# read counts vs GC content and variance?

indir <- "/Users/yeung/data/scchic/from_cluster/count_mat_cell_sums_from_bam_nonoverlapping"
infs <- list.files(indir, pattern = "*.csv", full.names = TRUE)
count.sum <- lapply(infs, fread) %>%
  bind_rows()  %>%
  rowwise() %>%
  mutate(mark = strsplit(cell, "-")[[1]][[4]]) 


# plog GC vs signal 
dat.impute.long <- melt(dat.impute.log)
colnames(dat.impute.long) <- c("coord", "cell", "log2exprs.imputed")

# add GC content per bin 
gr.gc.dat$coord <- paste(gr.gc.dat$chromo, paste(as.integer(gr.gc.dat$start), as.integer(gr.gc.dat$end), sep = "-"), sep = ":")

gr.gc.dat.merge <- left_join(gr.gc.dat, dat.impute.long)

gr.gc.dat.merge$experi <- sapply(as.character(gr.gc.dat.merge$cell), function(x) ClipLast(x, jsep = "_"))

gr.gc.dat.merge.sum <- gr.gc.dat.merge %>%
  group_by(experi, coord) %>%
  summarise(log2exprs.imputed = mean(log2exprs.imputed),
            gc = unique(gc))

ggplot(gr.gc.dat.merge, aes(x = gc, y = log2exprs.imputed)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi, ncol = 1)

ggplot(subset(gr.gc.dat.merge.sum, gc > 0.35), aes(x = gc, y = log2exprs.imputed)) + geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi, ncol = 1)

# plot GC vs H3K9me3 slopes

jsub.fits.merge <- left_join(jsub.fits, gr.gc.dat)

ggplot(jsub.fits.merge %>% filter(pos > 4e6), aes(x = gc, y = slope)) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + geom_point(alpha = 0.2) + geom_smooth(method = "lm", se = FALSE)

# do a heatmap but scaled to integrate both gc and the overall marks

jsub.fits.scaled <- jsub.fits %>% group_by(mark) %>% 
  # filter(pos > pos.filt) %>%
  mutate(signal = scale(slope)) %>%
  dplyr::select(pos, mark, signal)
gr.gc.scaled <- gr.gc.dat %>% 
  # filter(pos > pos.filt) %>%
  mutate(gc = ifelse(gc < quantile(gc, 0.01), quantile(gc, 0.01), gc)) %>%
  mutate(mark = "GC", signal = scale(gc)) %>%
  dplyr::select(pos, mark, signal)

jsub.merged.scaled <- bind_rows(jsub.fits.scaled, gr.gc.scaled)

jsub.merged.scaled$mark <- factor(as.character(jsub.merged.scaled$mark), levels = c(jmarks, "GC"))
ggplot(jsub.merged.scaled, aes(x = pos / 10^6, y = reorder(mark, dplyr::desc(mark)), fill = signal)) + geom_tile() + 
    theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") +
    scale_fill_gradientn(colours = colorRamps::matlab.like(5), name = "Signal") +
    ylab("") +
    xlab("Position (MB)")

# zooom into GC rich region
pos.start <- 15e6
pos.end <- 25e6
ggplot(jsub.merged.scaled %>% filter(pos > pos.start & pos < pos.end), aes(x = pos / 10^6, y = reorder(mark, dplyr::desc(mark)), fill = signal)) + geom_tile() + 
    theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") +
    scale_fill_gradientn(colours = colorRamps::matlab.like(5), name = "Signal") +
    ylab("") +
    xlab("Position (MB)")
# do correlatioins within this region?
# ggpairs(jsub.wide %>% filter(pos > pos.start & pos < pos.end))

m.cor.filt <- ggpairs(jsub.wide %>% filter(pos > pos.start & pos < pos.end) %>% dplyr::select(-pos),
                 lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5))) + theme_classic() + ggtitle("Filt")
print(m.cor.filt)

# dev.off()
# Jake Yeung
# Date of Creation: 2019-04-22
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/tf_activity_correlation_gene_levels.R
# Gene levels correlation

rm(list=ls())

tstart <- Sys.time()

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(JFuncs)
library(umap)
library(ggrepel)
library(tidyr)
library(hash)
library(tidytext)
library(ggrepel)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/TrajFunctions.R")
source("scripts/Rfunctions/AnnotationFunctions.R")

inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"

# Functions ---------------------------------------------------------------



# Constants ---------------------------------------------------------------

vcutoff <- 2
use.peaks.mat <- TRUE

pcutoff <- -7.5

jscale.fac <- 10^6
jpseudo <- 0
jsize <- 1
# make.plots <- FALSE
make.plots <- FALSE



# Load LDA exprs matrix ---------------------------------------------------

# tssdist <- 80000; jdate <- "2019-04-20"; Kvec <- "50"
# tssdist <- 60000; jdate <- "2019-04-22"; Kvec <- "50"
tssdist <- 50000; jdate <- "2019-04-22"; Kvec <- "50"
# tssdist <- 20000; jdate <- "2019-04-21"; Kvec <- "50"
# tssdist <- 40000; jdate <- "2019-04-20"; Kvec <- "25_50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_peaks_exprs_merge_clusters_explore/lda_out_meanfilt.PZ-BM-H3K4me1.CountThres0.K-", Kvec, "_GeneTSS.Dedup.", jdate, ".", tssdist, ".Robj")

assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

out.lda <- out.lda[[length(out.lda)]]
out.lda@documents <- SwitchColnames(unname(out.lda@documents), jsplit = "-")
tm.result <- posterior(out.lda)

top.cells <- tidy(out.lda, matrix = "gamma") %>%
  group_by(topic) %>%
  arrange(desc(gamma)) %>%
  mutate(rnk = seq(length(gamma))) %>%
  mutate(gamma.zscore = scale(gamma, center = TRUE, scale = TRUE)) %>%
  dplyr::rename(cell = document)

top.cells.sum <- top.cells %>%
  group_by(topic) %>% # do entropy on 1 to 99% of cells
  filter(gamma.zscore < quantile(gamma.zscore, 0.95)) %>%
  mutate(zscore.prob = exp(gamma.zscore) / sum(exp(gamma.zscore))) %>%
  summarise(entropy = -sum(zscore.prob * log(zscore.prob))) %>%
  arrange(entropy)

top.peaks <- tidytext::tidy(out.lda, matrix = "beta", log = FALSE) %>%
  group_by(topic) %>%
  arrange(desc(beta)) %>%
  mutate(rnk = seq(length(beta))) %>%
  mutate(beta.zscore = scale(beta, center = TRUE, scale = TRUE)) %>%
  rowwise() %>%
  mutate(gene = strsplit(term, ";")[[1]][[2]])

mat.impute <- t(tm.result$topics %*% tm.result$terms)

# get gene list
rnames <- rownames(mat.impute)
rnames.keep <- grepl(";", rnames)

mat.impute.sub <- mat.impute[rnames.keep, ]

genes <- sapply(rownames(mat.impute.sub), function(x){
  g <- tryCatch({
    return(strsplit(x, ";")[[1]][[2]])
  }, error = function(e){
    return("Peak")
  })
}, USE.NAMES = FALSE)

exprs.long <- data.frame(peak = rownames(mat.impute.sub), gene = genes, as.data.frame(mat.impute.sub)) %>%
  tidyr::gather(key = "cell", value = "exprs", c(-peak, -gene))


# Load data  --------------------------------------------------------------

load(inf, v=T)
load(inf.trajs, v=T)


# load MARA 

experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
switch.rep.hash <- GetRepSwitchHash(experihash)


indir.mara <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50"

mara.out <- LoadMARA(indir.mara, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = "")


# Plot UMAP with activities -----------------------------------------------

jmark <- "H3K4me1"
dat.sub <- dat.trajs.long %>% filter(mark == jmark)

# add to act long
act.exprs.umap <- left_join(mara.out$act.long, dat.sub %>% dplyr::select(X1, X2, mark, cell, louvain))


system.time(
  act.exprs.umap <- left_join(act.exprs.umap, exprs.long %>% dplyr::rename(motif = gene))
)

act.exprs.umap$exprs.log <- log10(act.exprs.umap$exprs * jscale.fac + jpseudo)

colact <- c("gray85", "gray50", scales::muted("darkred"))
colexprs <- c("gray85", "gray50", scales::muted("darkblue"))
jmotif <- "Tal1"
m1 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jcol = colact[[3]], jtitle = paste(jmotif, "activity"))
m2 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "exprs", jcol = colexprs[[3]], jtitle = paste(jmotif, "exprs"))
multiplot(m1, m2, cols = 2)

jmotif <- "Ebf1"
m1 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jcol = colact[[3]], jtitle = paste(jmotif, "activity"))
m2 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "exprs", jcol = colexprs[[3]], jtitle = paste(jmotif, "exprs"))
multiplot(m1, m2, cols = 2)

jmotif <- "Cebpb"
m1 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jcol = colact[[3]], jtitle = paste(jmotif, "activity"))
m2 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "exprs", jcol = colexprs[[3]], jtitle = paste(jmotif, "exprs"))
multiplot(m1, m2, cols = 2)

# 
# plot correlations??
jmotif <- "Ebf1"
jmotif <- "Cebpe"
jmotif <- "Cebpd"
jmotif <- "Cebpg"
jmotif <- "Cebpb"

jmotif <- "Tal1"
jmotif <- "Cebpb"
jmotif <- "Foxc1"
jmotif <- "Spic"




jmotif <- ""
m3 <- ggplot(act.exprs.umap %>% filter(motif == jmotif), aes(y = exprs.log, x = activity)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(jmotif)  + 
  geom_hline(yintercept = vcutoff) 
print(m3)

# 
# # Genome wide correlations ------------------------------------------------
# 
# act.sum <- act.exprs.umap %>%
#   group_by(motif) %>%
#   summarise(cor.out = cor(x = exprs.log, y = activity)) %>%
#   arrange(desc(abs(cor.out))) %>%
#   left_join(mara.out$zscores) %>%
#   mutate(motif.lab = ifelse(zscore > 2, motif, NA))
# 
# ggplot(act.sum, aes(x = cor.out, y = zscore, label = motif.lab)) + geom_point() + geom_text_repel() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + ggtitle(paste("Global analysis of TF and gene correlation:", tssdist, "TSS window"))
# 
# 
# # Add cell trajectories  --------------------------------------------------

ctypes <- c("granu", "eryth", "lymphoid")
traj.annot <- lapply(ctypes, function(jtraj) data.frame(cell = trajs.spring[[jmark]][[jtraj]]$cell, traj = jtraj, lambda = trajs.spring[[jmark]][[jtraj]]$lambda)) %>%
  bind_rows()

act.exprs.umap.withtraj <- left_join(act.exprs.umap, traj.annot) %>%
  filter(!is.na(traj))

jmotif <- "Cebpb"
jmotif <- "Foxc1"

jmotifs <- mara.out$zscores$motif[1:25]
# .log <- TRUE

pdf(paste0("/tmp/mara_output.", tssdist, ".pdf"), useDingbats = FALSE)
for (jmotif in jmotifs){
  m1 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "activity", jcol = colact[[3]], jtitle = paste(jmotif, "activity"))
  
  m2.log <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "exprs.log", jcol = colexprs[[3]], jtitle = paste(jmotif, "exprs"))
  m2 <- PlotXYWithColor(act.exprs.umap %>% filter(motif == jmotif), xvar = "X1", yvar = "X2", cname = "exprs", jcol = colexprs[[3]], jtitle = paste(jmotif, "exprs"))
  
  m3.col <- ggplot(act.exprs.umap.withtraj %>% filter(motif == jmotif), aes(x = exprs.log, y = activity, color = traj, group = traj)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE, color = "black") + 
    geom_vline(xintercept = vcutoff)  + 
    ggtitle(jmotif)
  m3.col.merge <- ggplot(act.exprs.umap.withtraj %>% filter(motif == jmotif), aes(x = exprs.log, y = activity, color = traj)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE, mapping = aes(x = exprs.log, y = activity), inherit.aes = FALSE, color = "black") + 
    geom_vline(xintercept = vcutoff)  + 
    ggtitle(jmotif)
  multiplot(m1, m2, cols = 2)
  multiplot(m1, m2.log, cols = 2)
  print(m3.col)
  print(m3.col.merge)
}
dev.off()

# Genome wide correlations ------------------------------------------------

act.sum <- act.exprs.umap %>%
  group_by(motif) %>%
  summarise(cor.out = cor(x = exprs.log, y = activity)) %>%
  arrange(desc(abs(cor.out))) %>%
  left_join(mara.out$zscores) %>%
  mutate(motif.lab = ifelse(zscore > 2, motif, NA))

ggplot(act.sum, aes(x = cor.out, y = zscore, label = motif.lab)) + geom_point() + geom_text_repel() +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("Global analysis of TF and gene correlation:", tssdist, "TSS window"))


# Do analysis by trajectory -----------------------------------------------

# y axis: mean exprs, x axis: correlation
act.exprs.umap <- left_join(act.exprs.umap, traj.annot)

act.sum.bytraj <- act.exprs.umap %>%
  group_by(motif, traj) %>%
  filter(!is.na(exprs.log)) %>%
  summarise(exprs.upper = quantile(exprs.log, 0.95),
            exprs.mean = mean(exprs.log),
            act.mean = mean(activity),
            cor.out = cor(x = exprs.log, y = activity)) %>%
  filter(!is.na(traj)) %>% 
  filter(exprs.upper > vcutoff) %>%
  left_join(., mara.out$zscores) %>% 
  mutate(motif.lab = ifelse(zscore > 2, motif, NA))

ggplot(act.sum.bytraj, aes(x = cor.out, y = exprs.mean, label = motif.lab)) + geom_point(alpha=0.2) +  facet_wrap(~traj)  + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_text_repel() + ggtitle("TSS dist:", tssdist)

ggplot(act.sum.bytraj, aes(x = cor.out, y = act.mean, label = motif.lab)) + geom_point(alpha=0.2) +  facet_wrap(~traj)  + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_text_repel() + ggtitle("TSS dist:", tssdist)


# Some kind of global analysis of z-scores  -------------------------------




# # What is exprs of Cebpb? -------------------------------------------------
# 
# 
# dat <- fread("/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv", sep = "\t")
# colnames(dat) <- gsub(" ", "_", colnames(dat))
# 
# dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
#   group_by(Gene_ID) %>%
#   mutate(FPKM = replace_na(FPKM, 0)) %>%
#   mutate(logFPKM = log2(FPKM + 1),
#          zscore = scale(logFPKM, center = TRUE, scale = TRUE))
# 
# jgene <- "Cebpg"
# jgene <- "Cebpb"
# jgene <- "Cebpe"
# jgene <- "Cebpa"
# # jgenes <- c("Cebpa", "Cebpb", "Cebpc", "Cebpd", "Cebpe", "Cebpg")
# jgenes <- c("Tmem189", "A530013C23Rik", "Ube2v1")
# jgenes1 <- "Cebpd"
# jsub <- dat.long %>% filter(Gene_Name %in% jgenes1) %>% arrange(desc(logFPKM)) %>% mutate(CellType = factor(CellType, levels = CellType))
# jsub.levels <- levels(jsub$CellType)
# 
# jsub <- dat.long %>% filter(Gene_Name %in% jgenes) %>% arrange(desc(logFPKM)) %>% mutate(CellType = factor(CellType, levels = jsub.levels))
# 
# ggplot(jsub, aes(x = CellType, y = logFPKM)) + geom_bar(stat = "identity") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) + 
#   ggtitle("Cebpb region") + facet_wrap(~Gene_Name) 
# 
# 

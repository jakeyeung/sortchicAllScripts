# Jake Yeung
# Date of Creation: 2019-04-10
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/tf_activity_analysis_correlations.r
# description# 

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(JFuncs)
library(umap)

source("scripts/Rfunctions/MaraDownstream.R")


source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")



inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"


inf.pvals <- "/Users/yeung/data/scchic/robjs/TFactivity_zscore_analysis.RData"

# Load data  --------------------------------------------------------------

load(inf, v=T)

load(inf.pvals, v=T)

rm(count.mat.lst)

head(mara.outs$H3K4me1$zscore)


# Show top hits -----------------------------------------------------------

pvals.top <- pvals.long %>%
  group_by(mark) %>%
  dplyr::top_n(n = 10, wt = -log10pval) %>%
  arrange(log10pval)
print(dim(pvals.top))
print(split(pvals.top, f = pvals.top$mark))

print(data.frame(subset(pvals.long, mark == "H3K4me1") %>% arrange(log10pval)))


# Calculate correlations --------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

ref.mark <- "H3K4me1"
jmark <- ref.mark
jcolvec <- c("gray85", "gray50", "darkblue")

peaks.filt <- (annots.lst[[jmark]] %>% filter(distanceToTSS == 0))$region_coord
rows.i <- rownames(count.imputed.lst[[jmark]]) %in% peaks.filt
counts.mat.sub <- count.imputed.lst[[jmark]][rows.i, ]

jmotif <- "Cebpb"
jgene <- "Cebpb"

jmotif <- "Cebpb"
jgene <- "Cebpb"

# jmotif <- "Tal1"
# jgene <- "Tal1"

jmotif <- "Ebf1"
jgene <- "Ebf1"

jmotif <- "Foxc1"
jgene <- "Foxc1"

jmotif <- "Pax6"
jgene <- "Pax6"

jmotif <- "Hmbox1"
jgene <- "Hmbox1"

jmotif <- "Pou2f2"
jgene <- "Pou2f2"

out.sub <- GetPeaksFromGene(jgene, annots.lst[[ref.mark]], dist = 0)
(jpeak <- SelectBestPeak(out.sub$peaks, regions.annot, tm.result.lst[[ref.mark]])[[1]])
jscale.fac <- 1
jpseudo <- 10^-6
jsize <- 1

# merge activity and gene counts

act.sub <- dat.merged.lst[[jmark]] %>% filter(motif == jmotif) %>% select(cell, motif, activity)

exprs.sub <- count.imputed.lst[[jmark]][jpeak, ]

exprs.sub.dat <- data.frame(exprs = exprs.sub, cell = names(exprs.sub))

merged.sub <- left_join(act.sub, exprs.sub.dat)

merged.sub.filt <- merged.sub %>% filter(exprs > quantile(exprs, probs = 0.2))

ggplot(merged.sub, aes(x = exprs, y = activity)) + geom_point() + theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(paste(jmotif, jgene))

ggplot(merged.sub.filt, aes(x = exprs, y = activity)) + geom_point() + theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(paste(jmotif, jgene))

ggplot(merged.sub, aes(y = exprs, x = activity)) + geom_point() + theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
  ggtitle(paste(jmotif, jgene))



# what's the correlation?
cor(merged.sub.filt$activity, merged.sub.filt$exprs, method = "pearson")
cor(merged.sub.filt$exprs, merged.sub.filt$activity, method = "spearman")

m.peak <- PlotImputedPeaks2(tm.result.lst[[jmark]], jpeak, jmarks[[jmark]],  
                            # use.count.mat = count.mat.lst[[jmark]],
                            use.count.mat = NULL,
                            usettings=custom.settings.new.lst[[jmark]], 
                            gname = jgene,
                            jsize = jsize, jcolvec = jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo, y.axis.factor = -1)
m.motif <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]] %>% mutate(umap2 = -1 * umap2), mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec)
multiplot(m.peak, m.motif, cols = 2)

# correlate for all exprs
exprs.long <- data.frame(bin = rownames(counts.mat.sub), as.data.frame(counts.mat.sub), stringsAsFactors = FALSE) %>%
  tidyr::gather(key = "cell", value = "exprs", -bin) %>%
  left_join(., act.sub) %>%
  filter(!is.na(activity))

# add motif of interest

# check
qplot(subset(exprs.long, bin == jpeak)$activity, subset(exprs.long, bin == jpeak)$exprs)
 
exprs.cors <- exprs.long %>%
  group_by(bin) %>% 
  # filter(exprs > quantile(exprs, probs = 0.2)) %>%
  summarise(cor.pears = cor(activity, exprs)) %>%
  arrange(desc(abs(cor.pears)))
print(head(exprs.cors))
# print(head(exprs.cors))

ggplot(exprs.cors, aes(x = cor.pears)) + geom_histogram(bins = 100) + ggtitle(paste(jmotif, jgene)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlim(c(-1, 1))

# plot positive or negative correlations

# jbin <- as.character(subset(exprs.cors, cor.pears == max(cor.pears))$bin)
# jbin <- as.character(subset(exprs.cors, cor.pears == max(cor.pears))$bin)
(jbin <- as.character(subset(exprs.cors, cor.pears == min(cor.pears))$bin))
jbin <- jpeak
# get gene form bin
(jgene <- subset(annots.lst[[jmark]], region_coord == jbin)$SYMBOL)
m.peak <- PlotImputedPeaks2(tm.result.lst[[jmark]], jbin, jmarks[[jmark]],  
                            # use.count.mat = count.mat.lst[[jmark]],
                            use.count.mat = NULL,
                            usettings=custom.settings.new.lst[[jmark]], 
                            gname = jgene,
                            jsize = jsize, jcolvec = jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo, y.axis.factor = -1)
m.motif <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark]] %>% mutate(umap2 = -1 * umap2), mara.outs[[jmark]]$zscores, jmark, jsize = 0.75, colvec = jcolvec)
multiplot(m.peak, m.motif, cols = 2)


# Jake Yeung
# Date of Creation: 2019-04-04
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/tf_activity_analysis_with_background.R
# Background of TF activity 
# 

library(data.table)
library(dplyr)
library(ggplot2)

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient.withColnameList.2019-04-04.RData"


# Load data  --------------------------------------------------------------

load(inf, v=T)

head(mara.outs$H3K4me1$zscore)

# Load zscore mats --------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

inf.zscores <- lapply(jmarks, function(jmark) paste0("/Users/yeung/data/scchic/from_cluster/zscore_permute_summary/", jmark, ".zscore_permute_summary.txt.gz"))

dats <- lapply(inf.zscores, function(inf.zscore) read.table(gzfile(inf.zscore), header = FALSE, stringsAsFactors = FALSE))

zscore.long <- lapply(dats, function(dat){
  colnames(dat) <- c("motif", "zscore", "seed", "mark")
  return(dat)
})

zscore.long <- do.call(rbind, zscore.long)

print(head(zscore.long))
print(unique(zscore.long$motif))

zscore.long <- subset(zscore.long, mark != "exiting")
subset(zscore.long, grepl("Zscores", motif))

zscore.long$zscore <- as.numeric(zscore.long$zscore)

zscore.long.real <- lapply(jmarks, function(jmark){
  return(mara.outs[[jmark]]$zscore %>% mutate(mark = jmark))
}) %>% 
  bind_rows() %>%
  rename(zscore.real = zscore)

zscore.long <- left_join(zscore.long, zscore.long.real)



# Plot CEBPB zscore versus background model  ------------------------------

GetPvalZscore <- function(jsub, zscore.real, jprob = 0.9, show.plot = TRUE, return.pval.only = FALSE, jtitle = ""){
  # zscore.real <- subset(zscore.real.dat, motif == jmotif)$zscore
  if (is.null(zscore.real)){
    # assume it's inside the jsub
    assertthat::assert_that(length(unique(jsub$zscore.real)) == 1)
    zscore.real <- jsub$zscore.real[[1]]
  }
  jsub <- jsub %>%
    arrange(zscore)
  # jsub$zscore <- as.numeric(jsub$zscore)
  jsub$zscore.cumsum <- cumsum(jsub$zscore)
  jsub$indx <- seq(nrow(jsub))
  jsub$frac.less.than <- jsub$indx / nrow(jsub)
  jsub$frac.more.than <- 1 - jsub$frac.less.than
  jsub$log10.frac.more.than <- log10(jsub$frac.more.than)
  
  jsubsub <- jsub %>% filter(zscore > quantile(zscore, probs = jprob) & frac.more.than > 0)
  jfit <- lm(formula = log10.frac.more.than ~ zscore, data = jsubsub)
  
  log10pval <- predict(jfit, newdata = data.frame(zscore = zscore.real))
  
  xpred <- seq(min(jsubsub$zscore), max(zscore.real, jsubsub$zscore), length.out = 100)
  ypred <- predict(jfit, newdata = data.frame(zscore = xpred))
  pred.dat <- data.frame(log10.frac.more.than = ypred, zscore = xpred)
  
  if (show.plot){
    m <- ggplot(jsubsub, aes(x = zscore, y = log10.frac.more.than)) + 
      geom_point() + theme_bw() +
      geom_vline(xintercept = zscore.real, linetype = "dashed") + 
      expand_limits(y = ceiling(log10pval)) + 
      geom_line(mapping = aes(x = zscore, y = log10.frac.more.than), data = pred.dat) + ggtitle(jtitle)
    print(m)
  }
  if (!return.pval.only){
    return(list(real.dat = jsubsub, pred.dat = pred.dat, fit = jfit, log10pval = log10pval))
  } else {
    return(data.frame(log10pval = log10pval))
  }
}

# out <- GetPvalZscore(zscore.long %>% filter(mark == jmark & motif == jmotif), subset(mara.outs$H3K4me1$zscore, motif == jmotif)$zscore)

pvals.long <- zscore.long %>%
  group_by(mark, motif) %>%
  do(GetPvalZscore(., zscore.real = NULL, jprob = 0.9, show.plot = FALSE, return.pval.only = TRUE))

# plot top hits by mark
ggplot(pvals.long, aes(x = log10pval)) + geom_histogram() + facet_wrap(~mark, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# show top hits

pvals.top <- pvals.long %>%
  group_by(mark) %>%
  dplyr::top_n(n = 10, wt = -log10pval) %>%
  arrange(log10pval)
print(dim(pvals.top))
print(split(pvals.top, f = pvals.top$mark))


# Save objects ------------------------------------------------------------

save(pvals.long, zscore.long, file = "~/data/scchic/robjs/TFactivity_zscore_analysis.RData")


# 
# out <- GetPvalZscore(zscore.long %>% filter(mark == jmark & motif == jmotif), zscore.real = NULL)
# 
# jprob <- 0.9
# 
# jmark <- "H3K4me1"
# jmotif <- "Nfatc1"
# jmotif <- "Cebpb"
# jmotif <- "Zbtb16"
# jmotif <- "Tal1"
# 
# jsub <- zscore.long %>% filter(mark == jmark & motif == jmotif) %>%
#   arrange(zscore)
# 
# zscore.real <- subset(mara.outs$H3K4me1$zscore, motif == jmotif)$zscore
# 
# jsub$zscore <- as.numeric(jsub$zscore)
# 
# ggplot(jsub, aes(x = log10(zscore))) + geom_density() + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_vline(xintercept = zscore.real) + 
#   ggtitle(jmotif)
# 
# # reverse cumulative?
# jsub$zscore.cumsum <- cumsum(jsub$zscore)
# jsub$indx <- seq(nrow(jsub))
# jsub$frac.less.than <- jsub$indx / nrow(jsub)
# jsub$frac.more.than <- 1 - jsub$frac.less.than
# jsub$log10.frac.more.than <- log10(jsub$frac.more.than)
# 
# # ggplot(jsub, aes(x = zscore)) + geom_density()
# ggplot(jsub, aes(x = zscore, y = log10(frac.more.than))) + 
#   geom_point() + theme_bw() +geom_vline(xintercept = zscore.real)
# ggplot(jsub %>% filter(zscore > quantile(zscore, probs = jprob)), aes(x = zscore, y = log10(frac.more.than))) + 
#   geom_point() + theme_bw() +geom_vline(xintercept = zscore.real)
# 
# # fit linear model
# jsubsub <- jsub %>% filter(zscore > quantile(zscore, probs = jprob) & frac.more.than > 0)
# jfit <- lm(formula = log10.frac.more.than ~ zscore, data = jsubsub)
# 
# log10pval <- predict(jfit, newdata = data.frame(zscore = zscore.real))
# 
# xpred <- seq(min(jsubsub$zscore), max(zscore.real, jsubsub$zscore), length.out = 100)
# ypred <- predict(jfit, newdata = data.frame(zscore = xpred))
# pred.dat <- data.frame(log10.frac.more.than = ypred, zscore = xpred)
# 
# ggplot(jsub %>% filter(zscore > quantile(zscore, probs = jprob)), aes(x = zscore, y = log10.frac.more.than)) + 
#   geom_point() + theme_bw() +
#   geom_vline(xintercept = zscore.real, linetype = "dashed") + 
#   expand_limits(y = ceiling(log10pval)) + 
#   geom_line(mapping = aes(x = zscore, y = log10.frac.more.than), data = pred.dat)


# Jake Yeung
# Date of Creation: 2019-05-06
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/1-quality_control_bins_and_cells.R
# Do quality control on bins and cells

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

fquad <- function(x, intercept, slope, curve){
  return(intercept + slope * x + curve * x^2)
}

GetMarkB6 <- function(x){
  return(strsplit(x, "-")[[1]][[4]])
}

GetRepB6 <- function(x){
  return(paste0("rep", strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
}

GetCellB6 <- function(x, shift = 0){
  # shift by 1 to get cells to start at 1 rather than 0
  cell <- as.numeric(strsplit(x, "_")[[1]][[4]]) + shift
  return(paste0("cell", as.character(cell)))
}

GetCellB6.uniqname <- function(x, shift = 0){
  # shift by 1 to get cells to start at 1 rather than 0
  cell <- as.numeric(strsplit(x, "_")[[1]][[2]]) + shift
  return(paste0("cell", as.character(cell)))
}

ReadDinuc <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName) 
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc) %>%
    rowwise() %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkB6(samp),
           repl = GetRepB6(samp),
           cell = GetCellB6.uniqname(samp, shift = 0))
  return(dat.long)
}

ReadDinucSplit <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName, strand = V2)
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc, -strand) %>%
    rowwise() %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkB6(samp),
           repl = GetRepB6(samp),
           cell = GetCellB6(samp))
  return(dat.long)
}

ReadMat <- function(inf, as.sparse = TRUE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  startend = V2)
  # remove missing
  dat <- subset(dat, startend != "Missing")
  dat$coord <- paste(dat$chromo, dat$startend, sep = ":")
  dat$chromo <- NULL
  dat$startend <- NULL
  if (as.sparse){
    coords <- dat$coord
    dat$coord <- NULL
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
    rownames(dat) <- coords
    dat <- Matrix::Matrix(dat, sparse = TRUE)
  }
  return(dat)
}

ReadCellSum <- function(inf){
  dat <- fread(inf)
}

# Load TA counts ----------------------------------------------------------

indir <- "/Users/yeung/data/scchic/from_cluster/TA_count_frequencies_bugfixed_uniqname"
# indir <- "/Users/yeung/data/scchic/from_cluster/TA_count_frequencies_bystrand_B6"
# indir <- "/Users/yeung/data/scchic/from_cluster/TA_count_frequencies-Ensembl95"

# infs <- list.files(indir, pattern = "*.RZcounts_bugfixed.test.csv", full.names = TRUE)
infs <- list.files(indir, pattern = "*.RZcounts_bugfixed.csv", full.names = TRUE)
# infs <- list.files(indir, pattern = "*.RZcounts.csv", full.names = TRUE)
# infs <- list.files(indir, pattern = "*.RZcounts_splitRS.csv", full.names = TRUE)

# dat.long <- lapply(infs, ReadDinucSplit) %>%
#   bind_rows()

system.time(
  dat.long <- lapply(infs, ReadDinuc) %>%
    bind_rows()
)

# Plot fraction of AT across samples, for each mark -----------------------

dat.sum <- dat.long %>%
  filter(dinuc != "Missing") %>%
  group_by(samp, mark) %>% 
  mutate(frac = count / sum(count)) %>%
  arrange(desc(count))

# dat.sum <- dat.long %>%
#   filter(dinuc != "Missing") %>%
#   group_by(samp, mark, strand) %>% 
#   mutate(frac = count / sum(count)) %>%
#   arrange(desc(count))

# ggplot(dat.sum %>% filter(dinuc == "AT"), aes(x = frac, fill = strand)) + geom_density(alpha = 0.2) + facet_wrap(~mark) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.sum %>% filter(dinuc == "AT"), aes(x = frac)) + geom_density(alpha = 0.2) + facet_wrap(~mark) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.sum %>% filter(dinuc == "TA"), aes(x = frac)) + geom_density(alpha = 0.2) + facet_wrap(~mark) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load total counts -------------------------------------------------------

indir <- "/Users/yeung/data/scchic/from_cluster/count_mat_cell_sums_from_bam"
infs <- list.files(indir, pattern = "*.csv", full.names = TRUE)

dat.sizes <- lapply(infs, ReadCellSum) %>%
  bind_rows()  %>%
  dplyr::rename(samp = cell)


# Label empty wells -------------------------------------------------------

# #  1 index
# indx.all <- seq(384)
# hascell.indx <- c(seq(1:356),seq(360:379)+360)
# empty.indx <- setdiff(indx.all, hascell.indx)
# empty.names <- paste("cell", empty.indx, sep = "")
# print(empty.names)

# 0 index
indx.all <- seq(384) - 1
hascell.indx <- c(seq(1:356),seq(360:379)+360) - 1
empty.indx <- setdiff(indx.all, hascell.indx)
empty.names <- paste("cell", empty.indx, sep = "")


# Scatterplot -------------------------------------------------------------

dat.merge <- left_join(dat.sizes, dat.sum %>% filter(dinuc == "TA"))
# dat.merge <- left_join(dat.sizes, dat.sum %>% filter(dinuc == "GC"))

dat.merge <- dat.merge %>%
  mutate(is.empty = cell %in% empty.names,
         cellsum.log = log10(cellsum + 1))

ggplot(dat.merge, aes(x = cellsum, y = frac, color = is.empty, size = is.empty)) + geom_point(alpha = 0.2) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(dat.merge, aes(x = cellsum + 1, y = frac, color = is.empty, size = is.empty)) + geom_point(alpha = 0.2) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_log10() + 
  facet_wrap(~mark) + xlab("# Unique Cuts") + ylab("Fraction of cuts starting with TA") 
  #scale_x_continuous(breaks= scales::pretty_breaks())

dat.means <- dat.merge %>%
  group_by(is.empty, mark) %>%
  summarise(cellsum.med = median(cellsum),
            cellsum.med.log = log10(cellsum.med + 1))

ggplot(dat.merge, aes(x = cellsum.log, fill = is.empty)) + geom_histogram(alpha = 0.5, bins = 50) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + xlab("log10(Unique Cuts + 1)") + geom_vline(mapping = aes(xintercept = cellsum.med.log), data = dat.means, linetype = "dotted")


# Separate empty versus nonempty wells with sigmoid  ----------------------

dat.merge.filt <- subset(dat.merge, !is.nan(frac))
jweights <- sapply(dat.merge.filt$is.empty, function(x) ifelse(x == TRUE, 50, 1))
# jweights.norm <- jweights / sum(jweights)
# glm.fit <- glm(is.empty ~ cellsum.log + I(cellsum.log ^ 2) + frac, data = dat.merge.filt, family = binomial, weights = jweights.norm)
# glm.fit <- glm(is.empty ~ cellsum.log + frac, data = dat.merge.filt, family = binomial, weights = jweights)
glm.fit <- glm(is.empty ~ cellsum.log + I(cellsum.log ^ 2) + frac, data = dat.merge.filt, family = binomial, weights = jweights)
# glm.fit <- glm(is.empty ~ cellsum.log + frac, data = dat.merge.filt, family = binomial)

# jpred <- predict(glm.fit)
dat.merge.filt$pred.resp <- predict(glm.fit, newdata = dat.merge.filt, type = "response")
dat.merge.filt$pred <- as.logical(round(dat.merge.filt$pred.resp))

dat.merge.filt$correct <- dat.merge.filt$is.empty * dat.merge.filt$pred

ggplot(dat.merge.filt, aes(x = cellsum.log, y = frac, color = is.empty, size = is.empty)) + geom_point(alpha = 0.2) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(dat.merge.filt, aes(x = cellsum.log, y = frac, color = is.empty, size = pred)) + geom_point(alpha = 0.2) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

# https://stats.stackexchange.com/questions/6206/how-to-plot-decision-boundary-in-r-for-logistic-regression-model
(intercept <- -1 * coef(glm.fit)[["(Intercept)"]] / coef(glm.fit)[["frac"]])
(slope <- -1 * coef(glm.fit)[["cellsum.log"]] / coef(glm.fit)[["frac"]])
curve <- -1 * coef(glm.fit)[["I(cellsum.log^2)"]] / coef(glm.fit)[["frac"]]


ggplot(dat.merge, aes(x = cellsum.log, y = frac, color = is.empty)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = slope, intercept = intercept)

ggplot(dat.merge, aes(x = cellsum.log, y = frac, color = is.empty)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_function(fun=fquad, colour="red", args = list(intercept = intercept, slope = slope, curve = curve)) + 
  ylim(c(0, 1))


p <- 0.95
pfrac <- 0.95
dat.merge.filt <- dat.merge.filt %>%
  # group_by(is.empty) %>%
  # ungroup() %>%
  group_by(mark) %>%
  mutate(good.cellsize = ifelse(cellsum.log > quantile(cellsum.log, probs = (1-p)), TRUE, FALSE),
         good.frac = ifelse(frac > quantile(frac, probs = (1-pfrac)), TRUE, FALSE),
         good.cell = as.logical(good.cellsize * good.frac))

ggplot(dat.merge.filt, aes(x = cellsum.log, y = frac, color = good.cell, size = is.empty)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~mark) + 
  ylim(c(0, 1)) + 
  ggtitle(paste(p, pfrac))


good.cells <- dat.merge.filt %>%
  group_by(mark) %>%
  filter(good.cell) %>%
  summarise(n.cells = length(good.cell)) %>%
  mutate(type = "filt")

orig.cells <- dat.merge %>%
  group_by(mark) %>%
  summarise(n.cells = length(samp)) %>%
  mutate(type = "orig")

cells.summary <- left_join(good.cells, orig.cells, "mark") %>%
  mutate(frac = n.cells.x / n.cells.y)

print(cells.summary)

ggplot(dat.merge.filt, aes(x = cellsum.log, y = frac, color = is.empty, size = is.empty)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~mark) + 
  ylim(c(0, 1))


# Fit sigmoid on good cells -----------------------------------------------

dat.merge.filt$good.cell.bin <- as.numeric(dat.merge.filt$good.cell)

glm.fit <- glm(good.cell ~ cellsum.log + I(cellsum.log ^ 2) + frac, data = dat.merge.filt, family = binomial)

# jpred <- predict(glm.fit)
dat.merge.filt$pred.resp <- predict(glm.fit, newdata = dat.merge.filt, type = "response")
dat.merge.filt$pred <- as.logical(round(dat.merge.filt$pred.resp))

dat.merge.filt$correct <- dat.merge.filt$good.cell * dat.merge.filt$pred

# https://stats.stackexchange.com/questions/6206/how-to-plot-decision-boundary-in-r-for-logistic-regression-model
(intercept <- -1 * coef(glm.fit)[["(Intercept)"]] / coef(glm.fit)[["frac"]])
(slope <- -1 * coef(glm.fit)[["cellsum.log"]] / coef(glm.fit)[["frac"]])
curve <- -1 * coef(glm.fit)[["I(cellsum.log^2)"]] / coef(glm.fit)[["frac"]]


ggplot(dat.merge.filt, aes(x = cellsum.log, y = frac, color = good.cell)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  stat_function(fun=fquad, colour="red", args = list(intercept = intercept, slope = slope, curve = curve)) + 
  ylim(c(0, 1)) + facet_wrap(~mark)


# Write good cells to file  -----------------------------------------------

good.cells <- data.frame(cell = subset(dat.merge.filt, good.cell)$samp)
fwrite(good.cells, file = paste0("~/data/scchic/quality_control/good_cells.pcounts.", (1-p), ".pfrac.", (1-pfrac), ".txt"), col.names = TRUE)


# check cells with 100% frac
dat.merge.filt %>% arrange(cellsum) %>% filter(good.cell)
# subset(dat.merge.filt, frac == 1) %>% arrange(desc(cellsum))

ggplot(dat.merge, aes(y = cellsum.log, x = is.empty)) + geom_boxplot()
ggplot(dat.merge, aes(y = frac, x = is.empty)) + geom_boxplot()


# indir <- "/Users/yeung/data/scchic/from_cluster/count_mat_from_bam"
# 
# infs <- list.files(indir, pattern = "*slidewin.csv.gz", full.names = TRUE)
# 
# dat.mat <- lapply(infs, ReadMat)
# 
# dat.cellsums <- lapply(dat.mat %>% dplyr::select(-coord), function(x) colSums(x, na.rm = TRUE))



# Plot 2D -----------------------------------------------------------------




# 
# # Check fastq -------------------------------------------------------------
# 
# dat.fastq.tmp <- readRDS("/Users/yeung/data/scchic/robjs/dat.sum_parsed_barcode_diinuc_counts.rds")
# plot(density(subset(dat.fastq.tmp, dinuc == "TA")$dinuc.frac))
# 
# # Weird? Compare with fastq analysis --------------------------------------
# 
# inmain <- "/Users/yeung/data/scchic/from_cluster/count_summaries_B6/parsed2_filt_counts"
# infs <- list.files(inmain, pattern = "*.parsed.out.gz", full.names = TRUE)
# 
# dat.fastq <- lapply(infs, function(inf){
#   dat.tmp <- fread(inf, col.names = c("countsSampname", "barcode", "dinuc"))
#   # split first col
#   dat.tmp$counts <- as.numeric(sapply(dat.tmp$countsSampname, function(x) strsplit(x, " ")[[1]][[1]]))
#   dat.tmp$sampname <- sapply(dat.tmp$countsSampname, function(x) strsplit(x, " ")[[1]][[2]])
#   dat.tmp$countsSampname <- NULL
#   return(dat.tmp)
# }) %>%
#   bind_rows()
# 
# dat.fastq <- dat.fastq %>%
#   group_by(sampname) %>%
#   mutate(frac = counts / sum(counts),
#          mark = strsplit(sampname, "\\.")[[1]][[1]])
#  
# # plot fraction of TA
# ggplot(dat.fastq %>% filter(dinuc == "TA"), aes(x = frac)) + geom_density() + facet_wrap(~mark)
# 

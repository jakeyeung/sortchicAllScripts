# Jake Yeung
# Date of Creation: 2020-09-25
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/RP1_cells_ChICTAPS/1-explore_RP1_with_index.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/facs_index_data/RP1_FUCCI_2020-09-25/FACS_pseudotime_curve_plotted.with_hoescht.txt"
# outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/facs_index_data/RP1_FUCCI_2020-09-25/FACS_pseudotime_curve_plotted.with_hoescht.pdf"

# pdf(outpdf, useDingbats = FALSE)


hubpath <- "/home/jyeung/hub_oudenaarden"


indir <- file.path(hubpath, "jyeung/data/scChiC/facs_index_data/RP1_FUCCI_2020-09-25")

infs <- list.files(indir, pattern = "*.csv", full.names = TRUE)

inf <- infs[[1]]

dats.clean <- lapply(infs, function(inf){
  indx <- fread(inf)
  bname <- strsplit(basename(inf), ".csv")[[1]][[1]]
  
  colnames(indx)
  
  grep("561", colnames(indx), value = TRUE)
  cname.y <- "*[561] 585/29"
  cname.x <- "*[488] 530/40"
  
  cnames.keep <- grepl(pattern = "585/29|530/40|Well", colnames(indx))
  dat.clean <- indx[, ..cnames.keep]
  dat.clean$bname <- bname
  return(dat.clean)
}) %>%
  bind_rows()
# cnames.keep <- colnames(dats.clean)[grepl(pattern = "^\\*\\585/29|^\\*\\[530/40", colnames(dats.clean))]
cnames.keep <- grep(pattern = "^\\*\\[561|^\\*\\[488", colnames(dats.clean), value = TRUE)

assertthat::assert_that(length(cnames.keep) == 2)

print(cnames.keep)

dats.hoescht <- lapply(infs, function(inf){
  indx <- fread(inf)
  bname <- strsplit(basename(inf), ".csv")[[1]][[1]]
  
  colnames(indx)
  
  # grep("561", colnames(indx), value = TRUE)
  # cname.y <- "*[561] 585/29"
  # cname.x <- "*[488] 530/40"
  
  # cnames.keep <- grepl(pattern = "^\\[405\\]\\ 460\\/50\\ Area$|Well", colnames(indx))
  cnames.keep <- grepl(pattern = "405|Well", colnames(indx))
  dat.clean <- indx[, ..cnames.keep]
  dat.clean$bname <- bname
  return(dat.clean)
}) %>%
  bind_rows()

dats.clean <- left_join(dats.clean, dats.hoescht, by = c("Well", "bname"))

# ensure no dupes
dats.clean$bnameWell <- interaction(dats.clean$Well, dats.clean$bname, sep = ";")
dats.clean <- subset(dats.clean, !duplicated(bnameWell))

ggplot(dats.clean, aes(y = `*[561] 585/29`, x = `*[488] 530/40`)) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 


# Get chromo / spikein counts ---------------------------------------------

indir.chromo <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/merged_by_mark/countTablesAndRZr1only_ByChromo.NewFilters"
infs.chromo <- list.files(indir.chromo, pattern = "*NoChromo.again.csv", full.names = TRUE)

jchromos <- c(seq(19), "X", "Y")

jspike <- "J02459.1"
# jspike <- "J02459.1"
inf.chromo <- infs.chromo[[1]]


dats.spikein <- lapply(infs.chromo, function(inf.chromo){
  dat.chromo <- GetChromoCounts(inf.chromo, spikeinchromo = jspike, chromos.keep = jchromos) %>%
    filter(chromo == "1") %>%
    rowwise() %>%
    mutate(rowcoord = AddPlateCoordinates(samp)$rowcoord, 
           colcoord = AddPlateCoordinates(samp)$colcoord)
  dat.chromo$fname <- basename(inf.chromo)
  return(dat.chromo)
}) %>%
  bind_rows()

# Integrate with idnex data  ----------------------------------------------


WellToRowAndCol <- function(well, rtrn = c("row", "col", "both")){
  # A1 -> row1, column1
  myLetters <- toupper(letters[1:26])
  row.alphabet <- substr(well, start = 1, stop = 1)
  row.i <- match(x = row.alphabet, table = myLetters)
  col.i <- as.numeric(strsplit(well, split = row.alphabet)[[1]][[2]])
  if (rtrn == "both"){
    return(c(row.i, col.i))
  } else if (rtrn == "row"){
    return(row.i)
  } else if (rtrn == "col"){
    return(col.i)
  }
}

dats.clean2 <- dats.clean %>% 
  rowwise() %>%
  mutate(bname = gsub(" ", "_", bname),
         bname = gsub("K36MB", "K36_MB", bname), 
         markindx = strsplit(bname, split = "_")[[1]][[3]],
         platechr = strsplit(bname, split = "_")[[1]][[6]],
         plate = as.numeric(platechr), 
         rowcoord = WellToRowAndCol(Well, "row"),
         colcoord = WellToRowAndCol(Well, "col"))

print(unique(dats.clean2$markindx))

dats.spikein2 <- dats.spikein %>%
  rowwise() %>%
  mutate(mark = strsplit(experi, "-")[[1]][[6]],
         platechr = strsplit(experi, "-")[[1]][[7]],
         plateinit = as.numeric(gsub("^pl", "", platechr))) %>%
  group_by(mark) %>%
  mutate(plate = dense_rank(plateinit))

print(unique(dats.spikein2$mark))

library(hash)
marklst <- list("k27me3" = "K27", "k36me3" = "K36", "k9me3" = "K09")
markhash.fwd <- hash::hash(marklst)
markhash <- hash::invert(markhash.fwd)

dats.clean2 <- dats.clean2 %>%
  rowwise() %>%
  mutate(mark = AssignHash(markindx, markhash, null.fill = NA))

# integrate now
dats.merge <- left_join(dats.spikein2, dats.clean2, by = c("mark", "rowcoord", "colcoord", "plate"))

dats.merge.filt <- subset(dats.merge, spikeincounts > 0)

ggplot(dats.merge.filt, aes(y = `*[561] 585/29`, x = `*[488] 530/40`, color = log2(chromocounts / spikeincounts))) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  facet_wrap(~mark) +
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 

ggplot(dats.merge.filt, aes(y = `*[561] 585/29`, x = `*[488] 530/40`, color = log2(chromocounts))) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  facet_grid(plate~mark) +
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 

# remove bad cells

ggplot(dats.merge, aes(x = log10(chromocounts))) + 
  geom_density() +
  theme_bw() + 
  facet_grid(plate~mark) +
  geom_vline(xintercept = 3) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dats.merge.filt <- subset(dats.merge, chromocounts > 1000)

ggplot(dats.merge.filt, aes(x = log2(chromocounts / spikeincounts))) + 
  geom_density() +
  theme_bw() + 
  facet_grid(plate~mark) +
  geom_vline(xintercept = 3) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jmarks <- c("k27me3", "k36me3", "k9me3")
indx <- seq(4)

for (jmark in jmarks){
  for (i in indx){
    m <- ggplot(subset(dats.merge.filt, mark == jmark & plate == i), aes(y = `*[561] 585/29`, x = `*[488] 530/40`, color = log2(chromocounts))) + 
      geom_point(alpha = 0.2)  + 
      theme_bw() + 
      facet_grid(plate~mark) +
      scale_color_viridis_c() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_x_log10() + scale_y_log10() 
    print(m)
  }
}

m <- ggplot(dats.merge.filt, aes(y = `*[561] 585/29`, x = `*[488] 530/40`, 
                                                                     color = log2(chromocounts / spikeincounts))) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  facet_grid(plate~mark) +
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 
print(m)


# "" 


ggplot(dats.merge.filt, aes(x = spikeincounts, y = chromocounts)) + 
  geom_point() 
 
ggplot(dats.merge.filt, aes(x = spikeincounts / totalcounts)) + 
  geom_histogram() 


# Read from AvO  ----------------------------------------------------------

indir.spikeinsonly <- file.path(hubpath, "jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/spikeins_only/read_UMI_summary")
infs.spikeinsonly <- list.files(path = indir.spikeinsonly, pattern = "*.csv", full.names = TRUE)

x <- infs.spikeinsonly[[1]]
names(infs.spikeinsonly) <- infs.spikeinsonly

dats.spikeinsonly <- lapply(infs.spikeinsonly, function(x){
  jdat <- fread(x, header = FALSE)
  colnames(jdat) <- c("samp", "duped_spikeins", "spikeins")
  jdat.out <- data.frame(samp = jdat$samp, spikeinsonly = jdat$spikeins, stringsAsFactors = FALSE)
  return(jdat.out)
}) %>%
  bind_rows()

dats.merge2 <- as.data.frame(left_join(dats.merge, dats.spikeinsonly, by = "samp"))
# cells.keep <- unique(dat.merge2$samp)


ggplot(dats.merge2 %>% filter(chromocounts > 1000), aes(y = `*[561] 585/29`, x = `*[488] 530/40`, color = log2(chromocounts / spikeinsonly))) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  facet_grid(plate ~ mark) +
  scale_color_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 


ggplot(dats.merge2, aes(x = spikeinsonly, y = spikeincounts)) + 
  geom_point(alpha = 0.25)  + 
  # scale_x_log10() + 
  # scale_y_log10() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Spikein UMIs by mapping to only Lambda phage") + 
  ylab("Spikein UMIs by mapping to Human + Lambda phage genome")


ggplot(dats.merge2, aes(x = chromocounts / spikeinsonly, y = chromocounts / spikeincounts)) + 
  geom_point()  + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Spikein UMIs by mapping to only Lambda phage") + 
  ylab("Spikein UMIs by mapping to Human + Lambda phage genome")

# fit principal curve? 

ggplot(dats.merge2, aes(y = `*[561] 585/29`, x = `*[488] 530/40`)) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() + scale_y_log10() 

ggplot(dats.merge2, aes(y = log10(`*[561] 585/29`), x = log10(`*[488] 530/40`))) + 
  geom_point(alpha = 0.2)  + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# ggplot(dats.merge2, aes(y = log10(`[561] 585/29`), x = log10(`[488] 530/40`))) + 
#   geom_point(alpha = 0.2)  + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(jmat, aes(y = log10(`*[561] 585/29`), x = log10(`*[488] 530/40`))) + 
#   geom_point(alpha = 0.2)  + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


library(princurve)

jmat <- dats.merge2[, cnames.keep]
rownames(jmat) <- dats.merge2$samp
jmat <- as.matrix(log10(jmat[complete.cases(jmat), ]))
colnames(jmat) <- c("x1", "x2")
jdat <- as.data.frame(jmat)

pout <- principal_curve(jmat)

pout.dat <- data.frame(x1 = pout$s[, 1], x2 = pout$s[, 2], lambda.raw = pout$lambda, stringsAsFactors = FALSE) %>%
  ungroup() %>%
  mutate(lambda = lambda.raw / sum(lambda.raw))

pout.pred <- project_to_curve(x = jmat, s = pout$s, stretch = 0)

plot(pout.pred$s)

plot(pout$s)

# jdat$samp <- dats.merge2$samp
jdat$samp <- rownames(jdat)
jdat$lambda.raw <- pout$lambda
jdat$lambda <- jdat$lambda.raw / sum(jdat$lambda.raw)


# plot(density(log10(jdat$dist_ind)))
# jdat.filt <- subset(jdat, log10(dist_ind) < -5)

jperiod <- max(pout.dat$lambda)
w <- 2 * pi / jperiod
pout.dat$lambda.periodic <- w * pout.dat$lambda

library(PhaseHSV)

pout.dat$col <- sapply(pout.dat$lambda, function(x) hsv(h = PhaseToHsv(phase.vec = x, min.phase = 0, max.phase = jperiod), s = 1, v = 1))
jdat$col <- sapply(jdat$lambda, function(x) hsv(h = PhaseToHsv(phase.vec = x, min.phase = 0, max.phase = jperiod), s = 1, v = 1))

# plot(pout$s)
  
ggplot(mapping = aes(x = x1, y = x2)) + 
  geom_line(data = pout.dat, mapping = aes(color = lambda.raw), size = 2) + 
  scale_color_viridis_c() + 
  geom_point(data = as.data.frame(jdat), alpha = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(mapping = aes(x = x1, y = x2, color = lambda)) + 
  geom_point(data = jdat, alpha = 0.1) + 
  scale_color_viridis_c() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(mapping = aes(x = x1, y = x2)) + 
  geom_point(data = pout.dat, mapping = aes(color = col), size = 2) + 
  scale_color_identity() + 
  geom_point(data = as.data.frame(jmat), alpha = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(mapping = aes(x = x1, y = x2, color = col)) + 
  geom_point(data = jdat, size = 2, alpha = 0.25) + 
  geom_path(data = pout.dat %>% arrange(desc(x1)), size = 2, color = "grey85") + 
  scale_color_identity() + 
  # geom_point(data = as.data.frame(jmat), alpha = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Green") + ylab ("Red")

dats.merge2 <- left_join(dats.merge2, jdat, by = "samp")

# write to output

# add bname and Well

# make better cnames
dats.merge2.copy <- dats.merge2
colnames(dats.merge2.copy) <- make.names(colnames(dats.merge2.copy))
jdat.output <- left_join(jdat, dats.merge2.copy %>% dplyr::select(c(-x1, -x2, -lambda.raw, -lambda, -col)), by = "samp")
# left_join(., dats.hoescht, by = c("Well", "bname"))
# jdat.output <- left_join(jdat, dats.merge2 %>% dplyr::select(c(samp, bname, Well)), by = "samp") %>%
#   left_join(., dats.hoescht, by = c("Well", "bname"))

fwrite(jdat.output, file = outf, sep = "\t")
dev.off()


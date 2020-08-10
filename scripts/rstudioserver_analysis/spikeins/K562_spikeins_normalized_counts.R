# Jake Yeung
# Date of Creation: 2020-07-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_spikeins_normalized_counts.R
# Look at normalized counts

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

outrds <- paste0("/home/jyeung/data/from_rstudioserver/spikein_fits/spikein_fits_bins.", Sys.Date(), ".logncells.rds")
outpdf <- paste0("/home/jyeung/data/from_rstudioserver/spikein_fits/spikein_fits_bins.", Sys.Date(), ".logncells.pdf")

# Functions ---------------------------------------------------------------


# Load genomeiwde summaries -----------------------------------------------

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.2020-07-22.rds"

jsub.sum <- readRDS(inf.meta)
dat.meta <- subset(jsub.sum, select = c(samp, conc, spikeincounts, chromocounts, spikeinconc, ncells))


# Load different spikeins measure counts in each chromosome ------------------------------------------------

spikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

hubpath <- "/home/jyeung/hub_oudenaarden"

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo"

infs.lst <- list.files(indir, pattern = "*.csv", full.names = TRUE)
names(infs.lst) <- sapply(infs.lst, basename)


jmark <- "H3K4me3"
jconc <- "37U"


# jmarks.lst <- c("H3K4me3", "H3K27me3"); names(jmarks.lst) <- jmarks.lst
# jconcs.lst <- c("37U", "75U"); names(jconcs.lst) <- jconcs.lst

spikeinconc.vec <- sort(unique(jsub.sum$spikeinconc)); names(spikeinconc.vec) <- spikeinconc.vec
prep.vec <- c("37U", "75U"); names(prep.vec) <- prep.vec
jmarks.vec <- c("H3K4me3", "H3K27me3"); names(jmarks.vec) <- jmarks.vec

inf <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.countTable.ByChromo.csv"))

assertthat::assert_that(file.exists(inf))

# Load data  --------------------------------------------------------------

# dat.filt.chromo.long <- lapply(jmarks.lst, function(jmark){
#   print(jmark)
#   lapply(jconcs.lst, function(jconc){
#     print(jconc)
#     inf <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.countTable.ByChromo.csv"))
#     print(inf)
#     dat.filt.long <- GetChromoCounts(inf)
#     dat.filt.long$mark <- jmark
#     dat.filt.long$conc <- jconc
#     return(dat.filt.long)
#   })  %>%
#     bind_rows()
# }) %>%
#   bind_rows()
# 

mat.bins.lst.lst <- lapply(jmarks.vec, function(jmark){
  print(jmark)
  lapply(prep.vec, function(jconc){
    print(jconc)
    inf <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.countTable.csv"))
    assertthat::assert_that(file.exists(inf))
    mat <- ReadMatSlideWinFormat(inf, as.sparse = TRUE, sort.rnames = FALSE, add.chromo = TRUE)
    return(mat)
  })
}) 

# filter good rows --------------------------------------------------------

mincounts <- -1

pdf(outpdf, useDingbats = FALSE)
mat.bins.filt.lst.lst <- lapply(jmarks.vec, function(jmark){
  print(jmark)
  lapply(prep.vec, function(jprep){
    print(jprep)
    jmat <- mat.bins.lst.lst[[jmark]][[jprep]]
    
    plot(density(log2(rowMeans(as.matrix(jmat)))))
    abline(v = mincounts)
    
    rows.keep <- which(log2(rowMeans(jmat)) > mincounts)
    print(dim(jmat))
    jmat.filt <- jmat[rows.keep, ]
    print(dim(jmat.filt))
    return(jmat.filt)
  })
})
dev.off()

# Filter good cells  ------------------------------------------------------


# first for columns are bad
# jmark <- "H3K4me3"
# jconc <- "37U"
# gene.input <- mat.bins.lst.lst[[jmark]][[jconc]]

make.test <- FALSE

print("Beginning fits...")
system.time(
  fits.out <- lapply(jmarks.vec, function(jmark){
    print(jmark)
    lapply(prep.vec, function(jprep){
      print(jprep)
      lapply(spikeinconc.vec, function(jspikeinconc){
        print(jspikeinconc)
        jcells <- subset(jsub.sum, spikeinconc == jspikeinconc & conc == jprep & mark == jmark & ncells > 0)$samp
        mat.filt <- mat.bins.filt.lst.lst[[jmark]][[jprep]][, jcells]
        # make long and fit by rows
        jrows <- rownames(mat.filt); names(jrows) <- jrows
        # jrows <- jrows[1:100]
        
        if (make.test){
          
          jrow <- jrows[[888]]
          input.dat <- data.frame(cell = colnames(mat.filt), genecounts = mat.filt[jrow, ], stringsAsFactors = FALSE) %>%
            left_join(., dat.meta, by = c("cell" = "samp")) %>%
            rowwise() %>%
            mutate(ncells = log(ncells))
          
          jfit <- FitGlmRowSpikeins(input.dat, return.fit.obj = TRUE)
          jfit.dat <- FitGlmRowSpikeins(input.dat, return.fit.obj = FALSE)
          jpred <- data.frame(ypred.unadj = predict(jfit, input.dat, se.fit = FALSE), ncells = input.dat$ncells, chromocounts = input.dat$chromocounts, spikeinconc = input.dat$spikeinconc, spikeincounts = input.dat$spikeincounts, stringsAsFactors = FALSE) %>%
            mutate(ypred = ypred.unadj - log(spikeincounts))
          
          ggplot(input.dat, aes(x = ncells, y = log(genecounts/ spikeincounts))) + geom_point() + 
            theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_line(data = jpred, mapping = aes(x = ncells, y = ypred)) + 
            xlab("log(ncells)") + 
            ggtitle(paste(jmark, jprep, jspikeinconc, jrow))
          
          ggplot(input.dat, aes(x = ncells, y = genecounts/ spikeincounts)) + geom_point() + 
            theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_smooth(method = "lm", se = FALSE) + 
            ggtitle(paste(jmark, jprep, jspikeinconc, jrow))
          
          ggplot(input.dat, aes(x = ncells, y = genecounts/ chromocounts)) + geom_point() + 
            theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_smooth(method = "lm", se = FALSE) + 
            ggtitle(paste(jmark, jprep, jspikeinconc, jrow))
          
          jfit <- FitGlmRowChromocounts(input.dat, return.fit.obj = TRUE)
          jfit.dat <- FitGlmRowChromocounts(input.dat, return.fit.obj = FALSE)
          jpred <- data.frame(ypred.unadj = predict(jfit, input.dat, se.fit = FALSE), ncells = input.dat$ncells, chromocounts = input.dat$chromocounts, spikeinconc = input.dat$spikeinconc, spikeincounts = input.dat$spikeincounts, stringsAsFactors = FALSE) %>%
            mutate(ypred = ypred.unadj - log(chromocounts))
          ggplot(input.dat, aes(x = ncells, y = log(genecounts / chromocounts))) + geom_point() + 
            theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_line(data = jpred, mapping = aes(x = ncells, y = ypred)) + 
            xlab("log(ncells)") + 
            ggtitle(paste(jmark, jprep, jspikeinconc, jrow))
          
        }
        
        
        
        # fit row
        lapply(jrows, function(jrow){
          input.dat <- data.frame(cell = colnames(mat.filt), genecounts = mat.filt[jrow, ], stringsAsFactors = FALSE) %>%
            left_join(., dat.meta, by = c("cell" = "samp")) %>%
            rowwise() %>%
            mutate(ncells = log(ncells))
          fit.dat <- FitGlmRowSpikeins(input.dat, return.fit.obj = FALSE)
          fit.dat$mark <- jmark
          fit.dat$prep <- jprep
          fit.dat$spikeinconc <- jspikeinconc
          return(fit.dat)
        })
      })
    })
  })
  
)
print("Beginning fits... done")

# Write fits to output ----------------------------------------------------


saveRDS(fits.out, file = outrds)





# Jake Yeung
# Date of Creation: 2020-07-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/gene_by_gene_fits_downstream.R
# Compare chromonorm vs spikein norm

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(ggrastr)


# Load genomeiwde summaries -----------------------------------------------

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.2020-07-22.rds"

jsub.sum <- readRDS(inf.meta)
dat.meta <- subset(jsub.sum, select = c(samp, conc, spikeincounts, chromocounts, spikeinconc, ncells)) %>%
  rowwise() %>%
  mutate(ncells.lin = ncells,
         ncells = log(ncells))




# Load fits -----------------------------------------------------------------

jdate <- "2020-07-24"
inf.spikein <- paste0("/home/jyeung/data/from_rstudioserver/spikein_fits/spikein_fits_bins.", jdate, ".logncells.FiltLow_800.lm.rds")
inf.chromo <- paste0("/home/jyeung/data/from_rstudioserver/spikein_fits/spikein_fits_bins.", jdate, ".NormByChromo.logncells.LowFilt_800.LM.rds")

assertthat::assert_that(file.exists(inf.spikein))
assertthat::assert_that(file.exists(inf.chromo))


pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/fits_gene_by_gene.", jdate, ".logncells.FiltLow.lm.pdf")


fits.spikein.lst <- readRDS(inf.spikein)
fits.chromo.lst <- readRDS(inf.chromo)

spikeinconc.vec <- sort(unique(jsub.sum$spikeinconc)); names(spikeinconc.vec) <- spikeinconc.vec
prep.vec <- c("37U", "75U"); names(prep.vec) <- prep.vec
jmarks.vec <- c("H3K4me3", "H3K27me3"); names(jmarks.vec) <- jmarks.vec


fits.spikein <- lapply(jmarks.vec, function(jmark){
  lapply(prep.vec, function(jprep){
    lapply(spikeinconc.vec, function(jconc){
      dat.tmp <- fits.spikein.lst[[jmark]][[jprep]][[as.character(jconc)]] %>%
        bind_rows()
    })  %>%
      bind_rows()
  }) %>%
    bind_rows()
})

fits.chromo <- lapply(jmarks.vec, function(jmark){
  lapply(prep.vec, function(jprep){
    lapply(spikeinconc.vec, function(jconc){
      dat.tmp <- fits.chromo.lst[[jmark]][[jprep]][[as.character(jconc)]] %>%
        bind_rows()
    })  %>%
      bind_rows()
  }) %>%
    bind_rows()
})



# jmark <- "H3K27me3"

pdf(pdfout, useDingbats = FALSE)

fits.lst <- list(spikeincounts = fits.spikein, chromocounts = fits.chromo)

jnames <- c("spikeincounts", "chromocounts")
xmax <- 4
ymax <- 20

jlab <- "Slope [log (counts) / log(cell)]"
for (jname in jnames){
  print(jname)
  subtitle <- paste("Norm By:", jname)
  
  # summarize all 
  m.sum <- ggplot(fits.lst[[jname]] %>% bind_rows() %>% mutate(spikeinfac = factor(spikeinconc, levels = unique(sort(spikeinconc)))), 
                  aes(x = spikeinfac, y = slope.ln)) + 
    geom_boxplot() + 
    facet_grid(mark ~ prep)  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "bottom") + 
    ylab(jlab) + 
    ggtitle(paste("Both marks. ", subtitle)) + 
    geom_hline(yintercept = c(0, 1), linetype = "dotted", color = 'blue') + 
    coord_cartesian(ylim = c(-xmax, xmax)) 
    # ylim(c(-xmax, xmax)) + 
  print(m.sum) 
  
  
  # compare slopes for K4me3 vs K27me3
  
  m4 <- ggplot(fits.lst[[jname]] %>% bind_rows(), aes(x = slope.ln, fill = mark)) + 
    geom_density(alpha = 0.25) + 
    facet_grid(prep ~ spikeinconc)  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.position = "bottom") + 
    # xlim(c(-xmax, xmax)) + 
    xlab(jlab) + 
    ggtitle(paste("Both marks. ", subtitle)) + 
    geom_vline(xintercept = c(0, 1), linetype = "dotted", color = 'blue') + 
    coord_cartesian(xlim = c(-xmax, xmax)) 
  print(m4) 


  for (jmark in jmarks.vec){
    print(jmark)
    # jfits.sub <- fits.spikein[[jmark]]
    jfits.sub <- fits.lst[[jname]][[jmark]]
    # (xmax <- quantile(abs(jfits.sub$slope), probs = 0.9999))
    
    m1 <- ggplot(jfits.sub, aes(x = slope, y = -log10(pval))) + 
      geom_point_rast(alpha = 0.2, size = 3) + 
      facet_grid(prep ~ spikeinconc)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      # xlim(c(-xmax, xmax)) + 
      xlab(jlab) + 
      geom_vline(xintercept = c(0, 1), linetype = "dotted", color = "blue") + 
      ggtitle(paste(jmark, subtitle)) + 
      coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax))  
    
    m2 <- ggplot(jfits.sub, aes(x = slope)) + 
      geom_density(alpha = 0.25, fill = "blue") + 
      facet_grid(prep ~ spikeinconc)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position = "bottom") + 
      geom_vline(xintercept = c(0, 1), linetype = "dotted", color = 'blue') + 
      ggtitle(paste(jmark, "distribution of slope estimates.", subtitle)) + 
      # xlim(c(-xmax, xmax)) + 
      coord_cartesian(xlim = c(-xmax, xmax)) + 
      xlab(jlab) 
    
    m3 <- ggplot(jfits.sub, aes(x = -log10(pval))) + 
      geom_density(alpha = 0.25, fill = "blue") + 
      facet_grid(prep ~ spikeinconc)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position = "bottom") + 
      ggtitle(paste(jmark, "distribution of pval estimates.", subtitle))
    
    m3.lin <- ggplot(jfits.sub, aes(x = pval)) + 
      geom_density(alpha = 0.25, fill = "blue") + 
      facet_grid(prep ~ spikeinconc)  + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         legend.position = "bottom") + 
      ggtitle(paste(jmark, "distribution of pval estimates linear.", subtitle)) 
    
    print(m1)
    print(m2)
    print(m3)
    print(m3.lin)
    
  }
}

# summarize for all conditions



dev.off()


# # plot theoretical
# 
# x <- seq(5)
# y <- x * 10
# 
# plot(log2(x), log2(y))
# 
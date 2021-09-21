# Jake Yeung
# Date of Creation: 2021-06-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/6-check_Wu_2021.R
# Recreate nreads boxplot

rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# PBMCs -------------------------------------------------------------------




inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/from_github/pbmc_fragments.tsv.gz"
dat.pbmc <- fread(inf)



colnames(dat.pbmc) <- c("chromo", "start", "end", "cell", "dupcount")
dat.pbmc$uniquecount <- 1

outrds.wu2021 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_total_reads_Wu2021_PBMCs.", Sys.Date(), ".rds")
saveRDS(dat.pbmc, file = outrds.wu2021k)

dat.pbmc.sum <- dat.pbmc %>%
  group_by(cell) %>%
  summarise(dupcount = sum(dupcount),
            uniquecount = sum(uniquecount))
  
ggplot(dat.pbmc.sum, aes(y = uniquecount)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.pbmc.sum, aes(x = dupcount, y = uniquecount)) + 
  geom_point() +
  scale_y_log10() + 
  scale_x_log10() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Check others ------------------------------------------------------------




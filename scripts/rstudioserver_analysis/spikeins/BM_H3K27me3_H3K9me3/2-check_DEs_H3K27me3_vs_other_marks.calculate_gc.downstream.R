# Jake Yeung
# Date of Creation: 2021-01-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/2-check_DEs_H3K27me3_vs_other_marks.calculate_gc.downstream.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(JFuncs)
library(scchicFuncs)



# Load it up --------------------------------------------------------------

load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_different_bin_sets.RData", v=T)

jnames <- names(gr.gc.dat.lst)
names(jnames) <- jnames

gr.gc.dat.long <- lapply(jnames, function(jname){
  jtmp <- gr.gc.dat.lst[[jname]]
  jtmp$jname <- jname
  return(jtmp)
}) %>%
  bind_rows()
  
ggplot(gr.gc.dat.long, aes(x = jname, y = gc)) + 
  geom_point() +
  scale_y_log10() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




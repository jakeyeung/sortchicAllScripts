# Jake Yeung
# Date of Creation: 2021-06-15
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/5-check_Ku_2021.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Load bed  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_total_reads_with_remapped_Ku2021.", Sys.Date(), ".AvgLibSizePerCell.pdf")

pdf(outpdf, useDingbats = FALSE)

# indir <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output")
# fnames <- list.files(indir, pattern = "*.bed.gz")
# fname.tmp <- "SRR10386929_1.bed.gz"
# # fname <- "SRR10615125_1.sorted_dupcounts.bed.gz"
# keeptop <- 30  # keep top 30 cells
# dat.sum <- lapply(fnames, function(fname){
#   print(fname)
#   inf <- file.path(indir, fname)
#   fbase <- strsplit(fname, split = "\\.")[[1]][[1]]
#   cnames <- c("chromo", "start", "end", "dupcounts", "bc")
#   dat.tmp <- fread(inf, header = FALSE, col.names = cnames)
#   dat.tmp$cellname <- paste(fbase, dat.tmp$bc, sep = "_") 
#   dat.tmp$uniquecounts <- 1
#   
#   dat.tmp.sum <- dat.tmp %>%
#     group_by(cellname) %>%
#     summarise(dupcounts = sum(dupcounts),
#               uniquecounts = sum(uniquecounts)) %>%
#     arrange(desc(uniquecounts))
#   # keep top 30 cellnames
#   dat.tmp.sum <- dat.tmp.sum[1:keeptop, ]
#   dat.tmp.sum$fbase <- fbase
#   return(dat.tmp.sum)
# })
# 
# dat.sum.long <- dat.sum %>%
#   bind_rows() 
# 

indir <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output/summarized_blacklist_filt")
infs <- list.files(path = indir, pattern = "*.bed", full.names = TRUE)
dat.sum <- lapply(infs, function(inf){
  fread(inf, header = TRUE)
})

dat.sum.long <- dat.sum %>%
  bind_rows()

keeptop <- 30
minreads <- 886
dat.sum.long <- lapply(dat.sum, function(jdat){
  rowskeep <- min(nrow(jdat), keeptop)
  jdat[1:rowskeep, ]
}) %>%
  bind_rows() %>%
  filter(uniquecounts > minreads)

ggplot(dat.sum.long, aes(y = uniquecounts / dupcounts)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.sum.long, aes(y = uniquecounts)) + 
  scale_y_log10() + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

outrds.ku2021 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_total_reads_Ku2021.", Sys.Date(), ".rds")
saveRDS(dat.sum.long, file = outrds.ku2021)

dat.sum.long$jset <- "Ku2021_remapped"

dat.toadd.uniq <- dat.sum.long %>%
  dplyr::select(cellname, uniquecounts, jset) %>%
  dplyr::rename(cell = cellname, 
                nreads = uniquecounts)

dat.toadd.uniqtotal <- dat.sum.long %>%
  dplyr::select(cellname, uniquecounts, jset, dupcounts) %>%
  dplyr::rename(cell = cellname, 
                nreads = uniquecounts,
                ntotal = dupcounts)

# Add to other data  ------------------------------------------------------

inf.rds <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.rds"
dat.other.uniq <- readRDS(inf.rds)

dat.merged.uniq <- rbind(dat.other.uniq, dat.toadd.uniq)

ggplot(dat.merged.uniq, aes(x = jset, y = nreads)) + 
  geom_boxplot() + 
  scale_y_log10() 

jcheck <- subset(dat.other.uniq, jset == "Ku_WBCs")

# Summaize  ---------------------------------------------------------------

inf.rds.uniqtotal <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds"
dat.other.uniqtotal <- readRDS(inf.rds.uniqtotal)
dat.merged.uniqtotal <- rbind(dat.other.uniqtotal, dat.toadd.uniqtotal)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.merged.uniqtotal, aes(x = ntotal, y = nreads, color = jset)) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(alpha = 0.25) +  
  scale_color_manual(values = cbPalette) + 
  scale_y_log10() + scale_x_log10()

ggplot(dat.merged.uniqtotal, aes(x = jset, y = nreads / ntotal)) + 
  theme_bw(24) + 
  ylab("Fraction of unique reads per sequenced read") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  geom_boxplot()  + 
  scale_color_manual(values = cbPalette) + 
  scale_y_log10()




dat.normfacs <- dat.merged.uniqtotal %>%
  group_by(jset) %>%
  summarise(meantotal = mean(ntotal),
            ncells = ifelse(unique(jset) != "Ku2021_remapped", length(cell), 9000),
            normfac = meantotal / ncells)

print(dat.normfacs)

dat.normfacs.add <- subset(dat.normfacs, jset == "Ku2021_remapped") %>%
  mutate(jset = "Ku_WBCs")
dat.normfacs.WithKu <- rbind(dat.normfacs, dat.normfacs.add)

dat.all.merge.fix.norm.WithKu <- left_join(dat.merged.uniq, dat.normfacs.WithKu)


# ggplot(dat.tmp.sum, aes(x = uniquecounts)) + 
#   geom_density()  + 
#   scale_x_log10() 

ggplot(dat.tmp.sum, aes(x = dupcounts, y = uniquecounts)) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()

ggplot(dat.all.merge.fix.norm.WithKu, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads)) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab("Unique reads") + 
  scale_y_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.all.merge.fix.norm.WithKu, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads / normfac)) + 
  ylab("Unique reads / AvgLibrarySizePerCell") + 
  geom_boxplot() + 
  theme_bw() + 
  scale_y_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()


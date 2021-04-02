# Jake Yeung
# Date of Creation: 2021-02-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_deeper/3-run_GREAT_on_regions.compare_GC.compare_dists.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(JFuncs)



hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/features_of_H3K27me3_quirky_regions"
outpdf <- file.path(outdir, paste0("H3K27me3_lost_flipped_bins_features.GO.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)


# Load bins  --------------------------------------------------------------

bsize <- 50000
inf.lost <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/compare_DE_analysis_spikein_vs_total/bins_lost.H3K27me3.2.576.2021-02-04.txt")
inf.flipped <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/compare_DE_analysis_spikein_vs_total/bins_flipped.H3K27me3.2.576.2021-02-04.txt")

dat.lost <- fread(inf.lost)
dat.flipped <- fread(inf.flipped)

# Run GREAT  --------------------------------------------------------------

regions.lost <- dat.lost %>%
  dplyr::rename(seqnames = Chr,
                start = Start,
                end = End)

regions.flipped <- dat.flipped %>%
  dplyr::rename(seqnames = Chr,
                start = Start,
                end = End)

regions.lost.mid <- dat.lost %>%
  dplyr::mutate(Start = Start + (bsize / 2) - 1,
                End = End - (bsize / 2) + 1) %>%
  dplyr::rename(seqnames = Chr,
                start = Start,
                end = End)

regions.flipped.mid <- dat.flipped %>%
  dplyr::mutate(Start = Start + (bsize / 2) - 1,
                End = End - (bsize / 2) + 1) %>%
  dplyr::rename(seqnames = Chr,
                start = Start,
                end = End)


regions.range.lost <- makeGRangesFromDataFrame(as.data.frame(regions.lost.mid))
regions.range.flipped <- makeGRangesFromDataFrame(as.data.frame(regions.flipped.mid))


regions.annotated.lost <- as.data.frame(annotatePeak(regions.range.lost,
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                annoDb='org.Mm.eg.db'))

regions.annotated.flipped <- as.data.frame(annotatePeak(regions.range.flipped,
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                annoDb='org.Mm.eg.db'))


gr.in.lost <- regions.range.lost
gr.in.flipped <- regions.range.flipped

out.great.lost <- submitGreatJob(gr.in.lost, species="mm10", request_interval = 30)
out.great.flipped <- submitGreatJob(gr.in.flipped, species="mm10", request_interval = 30)


out.tb.lost <- getEnrichmentTables(out.great.lost, ontology=availableOntologies(out.great.lost))
out.tb.flipped <- getEnrichmentTables(out.great.flipped, ontology=availableOntologies(out.great.flipped))



# Check lost bins  --------------------------------------------------------

out.tb.lost.ordered <- lapply(out.tb.lost, function(x){
  x %>% arrange(Hyper_Adjp_BH)
})

out.tb.flipped.ordered <- lapply(out.tb.flipped, function(x){
  x %>% arrange(Hyper_Adjp_BH)
})


lapply(out.tb.lost.ordered, function(x) print(head(x, n = 15)[, seq(3)]))
lapply(out.tb.flipped.ordered, function(x) print(head(x, n = 15)[, seq(3)]))


print(head(out.tb.lost.ordered$`GO Molecular Function`)[1:5, 1:5])
print(head(out.tb.flipped.ordered$`GO Molecular Function`)[1:5, 1:5])

# 
# head(out.tb.lost$`GO Molecular Function` %>% arrange(Hyper_Adjp_BH), n = 50)[, seq(4)]
# head(out.tb.lost$`GO Molecular Function` %>% arrange(Binom_Adjp_BH), n = 50)[, seq(4)]
# 
# head(out.tb.lost$`GO Biological Process` %>% arrange(Hyper_Adjp_BH), n = 50)[, seq(4)]
# head(out.tb.lost$`GO Biological Process` %>% arrange(Binom_Adjp_BH), n = 50)[, seq(4)]


# Plot GO -----------------------------------------------------------------


library(forcats)

out.tb.lst <- list("lost" = out.tb.lost.ordered, "flipped" = out.tb.flipped.ordered)
jnames <- names(out.tb.lst); names(jnames) <- jnames

m.lst.lst <- lapply(jnames, function(jname){
  out.tb <- out.tb.lst[[jname]]
  jterms <- names(out.tb); names(jterms) <- jterms
  m.lst.tmp <- lapply(jterms, function(jterm){
    jsub.go <- out.tb.lost.ordered[[jterm]][1:25, ]
    # clip long names
    jsub.go$name <- substr(x = jsub.go$name, start = 1, stop = 50)
    m <- ggplot(jsub.go, aes(x = forcats::fct_reorder(.f = name, .x = Hyper_Adjp_BH, .fun = median, .desc = FALSE), y = -log10(Hyper_Adjp_BH))) + 
      geom_col() + 
      theme_bw(12) + 
      ggtitle(paste(jname, jterm)) + 
      xlab("") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5))
    print(m)
    return(m)
  })
  return(m.lst.tmp)
})




# Filt by dist ------------------------------------------------------------



regions.annotated.lost.filt <- subset(regions.annotated.lost, abs(distanceToTSS) < 50000)
regions.annotated.flipped.filt <- subset(regions.annotated.flipped, abs(distanceToTSS) < 50000)


# Load K27me3 high bins as comparison  ------------------------------------



indir.bins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"

dat.bins <- lapply(jmarks, function(jmark){
  fname <- paste0("High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt")
  inf.bins <- file.path(indir.bins, fname)
  fread(inf.bins)
})

dat.bins <- lapply(dat.bins, function(jdat){
  jdat <- jdat %>%
    rowwise() %>%
    mutate(startExtend = start + 1 - bsize / 2,
           endExtend = end - 1 + bsize / 2)
  jdat$region_coordExtend <- paste(jdat$seqnames, paste(as.character(jdat$startExtend), as.character(jdat$endExtend), sep = "-"), sep = ":")
  return(jdat)
})

dat.de.bins <- lapply(jmarks, function(jmark){
  fname <- paste0("DE_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.", jmark, ".2021-01-30.txt")
  inf.bins <- file.path(indir.bins, fname)
  fread(inf.bins)
})

dat.de.bins <- lapply(dat.de.bins,function(jdat){
  jdat <- jdat %>%
    rowwise() %>%
    mutate(startExtend = start + 1 - bsize / 2,
           endExtend = end - 1 + bsize / 2)
  jdat$region_coordExtend <- paste(jdat$seqnames, paste(jdat$startExtend, jdat$endExtend, sep = "-"), sep = ":")
  return(jdat)
})


# Plot distance to gnes ---------------------------------------------------

plot(density(log10(abs(regions.annotated.flipped$distanceToTSS))))
plot(density(log10(abs(regions.annotated.lost$distanceToTSS))))

regions.dist.merged <- rbind(subset(regions.annotated.lost, select = c(seqnames, start, end, distanceToTSS, SYMBOL)) %>% mutate(Label = "Lost"),
                             subset(regions.annotated.flipped, select = c(seqnames, start, end, distanceToTSS, SYMBOL)) %>% mutate(Label = "Flipped"),
                             subset(dat.bins$H3K27me3, select = c(seqnames, start, end, distanceToTSS, SYMBOL)) %>% mutate(Label = "High"))

ggplot(regions.dist.merged, aes(x = log10(abs(distanceToTSS)), fill = Label)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(regions.dist.merged, aes(x = Label, y = log10(abs(distanceToTSS)), fill = Label)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Calculate GC content  ---------------------------------------------------

# calculaet GCs
# GCs
load("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_gc_analysis/gcs_genomewide.RData", v=T)

gr.gc.dat.dedup <- gr.gc.dat[!duplicated(gr.gc.dat), ]

plot(density(log10(abs(regions.annotated.flipped$distanceToTSS))))
plot(density(log10(abs(regions.annotated.lost$distanceToTSS))))


regions.gc.merged <- rbind(regions.lost %>% mutate(Label = "Lost"),
                           regions.flipped %>% mutate(Label = "Flipped"), 
                           subset(dat.bins$H3K27me3, select = c(seqnames, start, end, region_coordExtend)) %>% 
                                    dplyr::rename(Name = region_coordExtend) %>%
                                    mutate(Label = "High")) %>%
  left_join(., gr.gc.dat.dedup, by = c("Name" = "bname"))
                           
                     
ggplot(regions.gc.merged, aes(x = Label, y = gc, fill = Label)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Calculate GC content for the four marks  --------------------------------


regions.gc.merged.fourmarks <- lapply(jmarks, function(jmark){
  jtmp <- dat.bins[[jmark]]
  jtmp$Label <- jmark
  return(jtmp)
}) %>%
  bind_rows() %>%
  left_join(., gr.gc.dat.dedup, by = c("region_coordExtend" = "bname"))


ggplot(regions.gc.merged.fourmarks, aes(x = Label, y = gc, fill = Label)) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("High Bins") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(regions.gc.merged.fourmarks, aes(x = Label, y = gc)) + 
  geom_boxplot() + 
  theme_bw() + 
  ggtitle("High Bins") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Load pseudobulks HSPC only  ---------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})


# load LDA 
outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})


# spikeins only
dat.pbulk.lst <- lapply(jmarks, function(jmark.check){
  print(jmark.check)
  
  dat.metas.spikeins <- subset(dat.metas[[jmark.check]], !is.na(spikein_cuts)) %>%
    mutate(cluster = ifelse(cluster == "HSPCs", cluster, "notHSPCs"))
  
  cnames.keep.lst <- split(x = dat.metas.spikeins$cell, f = dat.metas.spikeins$cluster)
  pbulks.lst <- SumAcrossClusters(outs.all.lst[[jmark.check]]$count.mat, cnames.keep.lst = cnames.keep.lst)
  
  mat.pbulk <- bind_rows(pbulks.lst) %>%
    as.data.frame()
  rownames(mat.pbulk) <- rownames(outs.all.lst[[jmark.check]]$count.mat)
  mat.pbulk <- sweep(mat.pbulk, MARGIN = 2, STATS = colSums(mat.pbulk), FUN = "/")
  
  dat.pbulk <- as.matrix(mat.pbulk) %>%
    melt()
  colnames(dat.pbulk) <- c("bin", "ctype", "cuts")
  dat.pbulk$mark <- jmark.check
  return(dat.pbulk)
})

# annotate bins by lost 
dat.pbulk.hspcs.lost.lst <- lapply(dat.pbulk.lst, function(jdat){
  jdat.sub <- subset(jdat, ctype == "HSPCs") %>%
    rowwise() %>%
    mutate(is.lost = bin %in% regions.lost$Name)
})

dat.pbulk.nothspcs.lost.lst <- lapply(dat.pbulk.lst, function(jdat){
  jdat.sub <- subset(jdat, ctype != "HSPCs") %>%
    rowwise() %>%
    mutate(is.lost = bin %in% regions.lost$Name)
})

m.lst.hspcs <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.pbulk.hspcs.lost.lst[[jmark]], aes(x = log10(cuts * 1000000 + 1), fill = is.lost)) + 
    geom_density(alpha = 0.25) + 
    ggtitle(jmark, "HSPCs") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

m.lst.nothspcs <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.pbulk.nothspcs.lost.lst[[jmark]] %>% filter(ctype != "HSPCs"), aes(x = log10(cuts * 1000000 + 1), fill = is.lost)) + 
    geom_density(alpha = 0.25) + 
    ggtitle(jmark, "notHSPCs") + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

print(m.lst.hspcs)
print(m.lst.nothspcs)

m.lst <- lapply(jmarks, function(jmark){
  jmerge <- left_join(dat.pbulk.hspcs.lost.lst[[jmark]], dat.pbulk.nothspcs.lost.lst[[jmark]], by = c("bin")) %>%
    arrange(is.lost.x) 
  m <- ggplot(jmerge, aes(x = log2(cuts.x * 1000000 + 1), y = log2(cuts.y * 1000000 + 1), color = is.lost.x)) + 
    geom_point(alpha = 0.25) + 
    xlab("HSPCs signal") + ylab("nonHSPCs signal") + 
    ggtitle(jmark) + 
    theme_bw() + 
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})
print(m.lst)


dev.off()




# Write tables outputs ----------------------------------------------------

regions.gc.merged.fourmarks <- lapply(jmarks, function(jmark){
  jtmp <- dat.bins[[jmark]]
  jtmp$Label <- jmark
  return(jtmp)
}) %>%
  bind_rows() %>%
  left_join(., gr.gc.dat.dedup, by = c("region_coordExtend" = "bname"))

jlabs <- c("Flipped", "Lost")

for (jlab in jlabs){
  regions.gc.merged.annot.tmp <- left_join(regions.gc.merged %>% filter(Label == jlab) %>% dplyr::select(c(Name, Label, gc)), dat.bins$H3K27me3, by = c("Name" = "region_coordExtend"))
  outtxttmp <- file.path(outdir, paste0("H3K27me3_lost_flipped_bins_features.annotations.", jlab, ".", Sys.Date(), ".txt"))
  fwrite(regions.gc.merged.annot.tmp, file = outtxttmp, sep = "\t")
}


# 
# x <- 2^seq(10)
# scale(log(x))
# scale(log2(x))
# 


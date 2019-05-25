
ComparePublicLog <- function(inf, thres = 0.995, lab = "MyLabel", clip.chic.only = FALSE){
  dat <- data.table::fread(inf, sep = "\t")
  
  datref <- dat[, 1]
  colnames(datref) <- c("chip")
  datcompare <- dat[, 2:ncol(dat)]
  
  datref$peakid <- seq(nrow(datref))
  datcompare$peakid <- seq(nrow(datcompare))
  
  # plot 1 vs all
  datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")
  
  dat.long <- dplyr::left_join(datcompare.long, datref)
  
  dat.long.filt <- dat.long %>%
    filter(!is.nan(chip)) %>%
    group_by(Sample) %>%
    mutate(compare = lab)
  
  
  if (!clip.chic.only){
    dat.long.filt <- dat.long.filt %>%
      group_by(Sample) %>%
      filter(chic <= quantile(chic, probs = jthres) & chip <= quantile(chip, probs = jthres) & chic >= quantile(chic, probs = 1 - jthres) & chip >= quantile(chip, probs = 1-jthres))
  } else {
    dat.long.filt <- dat.long.filt %>%
      group_by(Sample) %>%
      filter(chic <= quantile(chic, probs = jthres) & chic >= quantile(chic, probs = 1-0.9*jthres))
  }
  dat.cors <- dat.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chip, chic, method = "pearson"),
              corr.spear = cor(chip, chic, method = "spearman")) %>%
    mutate(compare = lab)
  # remoev outliers?
  
  # m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(dat.long.filt, aes(x = chip, y = chic)) + 
    geom_hex(bins = 25) +
    # geom_point(alpha = 0.01) + 
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~Sample)
  return(list(dat.long.filt = dat.long.filt, dat.cors = dat.cors, m = m))
}

ComparePublicLinear <- function(inf, thres = 0.995, lab = "MyLabel", pseudocount = 1, filter.min = FALSE){
  dat <- data.table::fread(inf, sep = "\t")
  
  datref <- dat[, 1]
  colnames(datref) <- c("chip")
  datcompare <- dat[, 2:ncol(dat)]
  
  datref$peakid <- seq(nrow(datref))
  datcompare$peakid <- seq(nrow(datcompare))
  
  # plot 1 vs all
  datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")
  
  dat.long <- dplyr::left_join(datcompare.long, datref)
  
  # remove NaNs
  
  if (!filter.min){
    dat.long.filt <- dat.long %>%
      filter(!is.nan(chip)) %>%
      group_by(Sample) %>%
      filter(log(chic) <= quantile(log(chic), probs = thres) & log(chic) >= quantile(log(chic), probs = 1 - thres)) %>%
      filter(log(chip) <= quantile(log(chip), probs = thres) & log(chip) >= quantile(log(chip), probs = 1 - thres)) %>%
      mutate(compare = lab)
  } else {
    dat.long.filt <- dat.long %>%
      filter(!is.nan(chip)) %>%
      group_by(Sample) %>%
      # filter(chic > min(chic) & chip > min(chip)) %>%
      filter(chic > 5e-2) %>%
      filter(log(chic) <= quantile(log(chic), probs = thres) & log(chic) >= quantile(log(chic), probs = 1 - thres)) %>%
      filter(log(chip) <= quantile(log(chip), probs = thres) & log(chip) >= quantile(log(chip), probs = 1 - thres)) %>%
      mutate(compare = lab)
    
  }
  
  # dat.long.filt <- subset(dat.long.filt, !is.nan(chip))
  
  dat.cors <- dat.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chip, chic, method = "pearson"),
              corr.spear = cor(chip, chic, method = "spearman"),
              corr.pears.log = cor(log(chip + pseudocount), log(chic + pseudocount), method = "pearson"),
              corr.spear.log = cor(log(chip + pseudocount), log(chic + pseudocount), method = "spearman")) %>%
    mutate(compare = lab)
  # remoev outliers?
  
  # m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(dat.long.filt, aes(x = log10(chip), y = log10(chic))) + 
    geom_hex(bins = 25) +
    # geom_point(alpha = 0.01) + 
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~Sample)
  return(list(dat.long.filt = dat.long.filt, dat.cors = dat.cors, m = m))
}

CompareAcrossMarks <- function(inf2, jmark, refmark, clstr, thres = 0.995, lab = "MyLabel", clip.top.only = FALSE){
  dat.multi <- data.table::fread(inf2, sep = "\t")
  # compare erythryroblasts
  refname <- paste0(refmark, "_cluster_", clstr, ".bw")
  datref.multi <- dat.multi[, ..refname]
  # check we got the right cluster??
  print(head(datref.multi))
  colnames(datref.multi) <- c("chic_ref")
  cols.compare <- grepl(jmark, colnames(dat.multi))
  datcompare.multi <- dat.multi[, ..cols.compare]
  
  datref.multi$peakid <- seq(nrow(datref.multi))
  datcompare.multi$peakid <- seq(nrow(datcompare.multi))
  
  datcompare.multi.long <- melt(datcompare.multi, id.vars = "peakid", variable.name = "Sample", value.name = "chic")
  
  dat.multi.long <- dplyr::left_join(datcompare.multi.long, datref.multi)
  
  # remove outliers
  if (!clip.top.only){
    dat.multi.long.filt <- dat.multi.long %>%
      group_by(Sample) %>%
      filter(!is.nan(chic) & !is.nan(chic_ref)) %>%  # otherwise fails 
      filter(chic <= quantile(chic, probs = thres) & chic >= quantile(chic, probs = 1 - thres)) %>%
      filter(chic_ref <= quantile(chic_ref, probs = thres) & chic_ref >= quantile(chic_ref, probs = 1 - thres)) %>%
      mutate(compare = lab)
  } else {
    dat.multi.long.filt <- dat.multi.long %>%
      group_by(Sample) %>%
      filter(!is.nan(chic) & !is.nan(chic_ref)) %>%  # otherwise fails 
      filter(chic <= quantile(chic, probs = thres) & chic_ref <= quantile(chic_ref, probs = thres)) %>%
      mutate(compare = lab)
  }
  
  m <- ggplot(dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~Sample) + theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  dat.cors <- dat.multi.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chic, chic_ref, method = "pearson"),
              corr.spear = cor(chic, chic_ref, method = "spearman"),
              corr.pears.log2 = cor(log2(chic), log2(chic_ref), method = "pearson"),
              corr.spear.log2 = cor(log2(chic), log2(chic_ref), method = "spearman")) %>%
    
    mutate(compare = lab)
  # print(dat.cors)
  return(list(dat.multi.long.filt = dat.multi.long.filt, dat.cors = dat.cors, m = m))
}

CompareAcrossMarks2 <- function(inf2, jmark, clstr, thres = 0.995, lab = "MyLabel", clip.top.only = FALSE, jmark.ref){
  dat.multi <- data.table::fread(inf2, sep = "\t")
  # compare erythryroblasts
  clstr.name <- paste0(jmark.ref, "_cluster_", clstr, ".bw")
  datref.multi <- dat.multi[, ..clstr.name]
  colnames(datref.multi) <- c("chic_ref")
  cols.compare <- grepl(jmark, colnames(dat.multi))
  datcompare.multi <- dat.multi[, ..cols.compare]
  
  datref.multi$peakid <- seq(nrow(datref.multi))
  datcompare.multi$peakid <- seq(nrow(datcompare.multi))
  
  datcompare.multi.long <- melt(datcompare.multi, id.vars = "peakid", variable.name = "Sample", value.name = "chic")
  
  dat.multi.long <- dplyr::left_join(datcompare.multi.long, datref.multi)
  
  # remove outliers
  if (!clip.top.only){
    dat.multi.long.filt <- dat.multi.long %>%
      group_by(Sample) %>%
      filter(!is.nan(chic) & !is.nan(chic_ref)) %>%  # otherwise fails 
      filter(chic <= quantile(chic, probs = thres) & chic >= quantile(chic, probs = 1 - thres)) %>%
      filter(chic_ref <= quantile(chic_ref, probs = thres) & chic_ref >= quantile(chic_ref, probs = 1 - thres)) %>%
      mutate(compare = lab)
  } else {
    dat.multi.long.filt <- dat.multi.long %>%
      group_by(Sample) %>%
      filter(!is.nan(chic) & !is.nan(chic_ref)) %>%  # otherwise fails 
      filter(chic <= quantile(chic, probs = thres) & chic_ref <= quantile(chic_ref, probs = thres)) %>%
      mutate(compare = lab)
  }
  
  m <- ggplot(dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) + 
    # scale_x_log10() + scale_y_log10() + facet_wrap(~Sample) + theme_bw() + 
    facet_wrap(~Sample) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  dat.cors <- dat.multi.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chic, chic_ref, method = "pearson"),
              corr.spear = cor(chic, chic_ref, method = "spearman"),
              corr.pears.log2 = cor(log2(chic), log2(chic_ref), method = "pearson"),
              corr.spear.log2 = cor(log2(chic), log2(chic_ref), method = "spearman")) %>%
    mutate(compare = lab)
  # print(dat.cors)
  return(list(dat.multi.long.filt = dat.multi.long.filt, dat.cors = dat.cors, m = m, dat.multi.long = dat.multi.long, datref.multi = datref.multi, datcompare.multi = datcompare.multi))
}
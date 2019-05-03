
GetPvalFOV <- function(fovs.permute, fovs.real = NULL, jprob = 0.9, show.plot = FALSE, return.pval.only = FALSE){
  jmark <- fovs.permute$mark[[1]]
  if (!is.null(fovs.real)){
    fov.real <- subset(fovs.real, mark == jmark)$fov
  } else {
    assertthat::assert_that(length(unique(fovs.permute$fov.real)) == 1)
    fov.real <- fovs.permute$fov.real[[1]]
  }
  
  fovs.sub <- fovs.permute %>%
    # filter(mark == jmark) %>%
    group_by(mark) %>%
    arrange(fov) %>%
    ungroup() %>%
    mutate(fov.cumsum = cumsum(fov),
           indx = seq(length(fov)),
           frac.less.than = indx / length(indx),
           frac.more.than = 1 - frac.less.than,
           log10.frac.more.than = log10(frac.more.than))
  
  
  fovs.subsub <- fovs.sub %>% filter(fov > quantile(fov, probs = jprob) & frac.more.than > 0)
  jfit <- lm(formula = log10.frac.more.than ~ fov, data = fovs.subsub)
  
  log10pval <- predict(jfit, newdata = data.frame(fov = fov.real))
  
  xpred <- seq(min(fovs.subsub$fov), max(fov.real, fovs.subsub$fov), length.out = 100)
  ypred <- predict(jfit, newdata = data.frame(fov = xpred))
  pred.dat <- data.frame(log10.frac.more.than = ypred, fov = xpred)
  
  if (show.plot){
    m <- ggplot(fovs.subsub, aes(x = fov, y = log10.frac.more.than)) + 
      geom_point() + theme_bw() +
      geom_vline(xintercept = fov.real, linetype = "dashed") + 
      expand_limits(y = ceiling(log10pval)) + 
      geom_line(mapping = aes(x = fov, y = log10.frac.more.than), data = pred.dat) + ggtitle(jmark)
    print(m)
  }
  if (!return.pval.only){
    return(list(real.dat = fovs.subsub, pred.dat = pred.dat, fit = jfit, log10pval = log10pval, fov.real = fov.real))
  } else {
    return(data.frame(log10pval = log10pval))
  }
}



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
    return(list(real.dat = jsubsub, pred.dat = pred.dat, fit = jfit, log10pval = log10pval, zscore.real = zscore.real))
  } else {
    return(data.frame(log10pval = log10pval))
  }
}

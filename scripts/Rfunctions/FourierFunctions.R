# Jake Yeung
# Nov 10 2014
# functions related to Mr. Fourier


NGenesByAmp.long <- function(amps.long, amp.thres, labelid="tissue", varid="amp", outlabel = "tissue"){
  tiss <- amps.long[[labelid]][[1]]
  return(NGenesByAmp(amps.long[[varid]], amp.thres, tiss, outlabel))
}

NGenesByAmp <- function(amps, amp.thres, lab="", outlabel="tissue"){
  n.genes <- rep(length(amp.thres), times = length(amp.thres)) 
  for (i in seq(amp.thres)){
    amp.t <- amp.thres[i]
    n.genes[i] <- length(which(amps > amp.t))
  }
  out.dat <- data.frame(amp.thres = amp.thres, n.genes = n.genes)
  out.dat[[outlabel]] <- lab
  return(out.dat)
}

FracVarByGeneList <- function(dat.var.s, dat.var.filt, dat.complex.all_T, genes.sub, period.for.sort = 24, jlab = NULL){
  total.var.24.12 <- hash(paste(as.character(subset(dat.var.s, period %in% c(12, 24))$tissue), subset(dat.var.s, period %in% c(12, 24))$period, sep = ";"), subset(dat.var.s, period %in% c(12, 24))$sum_sqr_mod)
  dat.var.filt.sub <- subset(dat.var.filt, gene %in% genes.sub)
  dat.complex.all_T.sub <- subset(dat.complex.all_T, gene %in% genes.sub & period %in% c(24, 12)) %>%
    group_by(tissue, period) %>%
    summarise(var.sub = sum(Mod(exprs.transformed) ^ 2) * 2)
  dat.complex.all_T.sub$var.temp.total <- sapply(paste(as.character(dat.complex.all_T.sub$tissue), dat.complex.all_T.sub$period, sep = ";"), function(tiss) total.var.24.12[[tiss]])
  dat.complex.all_T.sub$var.norm <- dat.complex.all_T.sub$var.sub / dat.complex.all_T.sub$var.temp.total
  if (is.numeric(period.for.sort)){
    dat.sub <- subset(dat.complex.all_T.sub, period == period.for.sort)
    dat.complex.all_T.sub$tissue <- factor(dat.complex.all_T.sub$tissue, levels = dat.sub$tissue[order(dat.sub$var.norm, decreasing = TRUE)])
  }
  # add label
  if (!is.null(jlab)) dat.complex.all_T.sub$label <- jlab
  return(dat.complex.all_T.sub)
}


CalculatePeriodogram <- function(x, is.matrix=FALSE){
  # Creates periodogram for spectral analysis.
  # 
  # INPUT:
  # x = vector of data you want to check for rhythmicity
  # 
  # OUTPUT:
  # list with $freq being frequencies and periodogram values $periodogram
  #
  # Uses fast fourier transform
  # Adopted from https://onlinecourses.science.psu.edu/stat510/node/71
  
  # BEGIN: do my FFT, extract only relevant ranges
  if (is.matrix==FALSE){
    N <- length(x)  
  }
  else{
    N <- ncol(x)
  }
  # Create my "unscaled" periodogram
  
  if (is.matrix==FALSE){
    FF <- abs(fft(x) / sqrt(N)) ^ 2  # do my FFT. Fast Fourier Transform
  }
  else {
    FF <- t(mvfft(t(x)))
    # FF <- t(abs(mvfft(t(x))))
    # FF <- t(abs(mvfft(t(x)) / sqrt(N)) ^ 2)  # do my FFT. Fast Fourier Transform
  }
  # Create my "scaled" periodogram. Scale constant is 1/n
  scale.factor <- 1 / N
  # only need first (N/2) + 1 values of FFT result for periodogram
  
  if (is.matrix==FALSE){
    P <- scale.factor * FF[1:(N / 2 + 1)]
    P.unscaled <- FF[1:(N / 2 + 1)]    
  }
  else {
    P <- scale.factor * FF[, 1:(N / 2 )]
    P.unscaled <- FF[, 1:(N / 2 )]  
  }
  
  
  # creates harmonic frequencies from 0 to 0.5 in steps of 1/N for periodogram.
  # I only need from 0 to 0.5.
  f <- (0:(N / 2) / N)  
  # END: do my FFT, extract only relevant ranges
  
  # Why can't R return two objects? Returning list instead...
  return(list("freq"=f, "p.scaled"=P, "p.unscaled"=P.unscaled))
}

FindMaxFreqs <- function(freq, periodogram, n=5) {
  # Given periodogram and frequency, return frequency at which
  # maximum value of periodogram occurs
  # 
  # Input:
  # f = frequency calculated from CalculatePeriodogram
  # P = periodogram calculated from CalculatePeriodogram
  # n = return the top n frequencies. Default 5
  # 
  # Output:
  # max frequency
  # 
  max.vals <- sort(periodogram, decreasing=TRUE)[1:n]
  max.indices <- match(max.vals, periodogram)
  # Get freqs from indices
  max.freqs <- freq[max.indices]
  
  return(max.freqs)
}

PlotPeriodogram <- function(Frequency, Periodogram, title="Plot title", vline=NA, cex=1) {
  # Plots periodogram. 
  # 
  # Input: 
  # f = "frequency" calculated from CalculatePeriodogram
  # P = "periodogram" calculated from CalculatePeriodogram
  # 
  # Output:
  # Periodogram plot
  plot(Frequency, Periodogram, type='l', main=title, 
       cex.axis = cex,
       cex.lab = cex,
       cex.main = cex)
}

ProjectToPeriodicTime <- function(Y , N.TISSUES, N.TIMEPTS, INTERVAL, OMEGA, col.names){
  # Project matrix of times to frequency domain
  # ARGS:
  #   Y: matrix of expression. Rows are genes. Columns contain tissues and time.
  #      In this case, expect tissues clustered together, ordered by time.
  #   col.names: column names of output Y's projected onto time. FALSE means no col.names
  #   N.TISSUES: number of tissues
  #   N.TIMEPTS: number of time points per tissue
  #   INTERVAL: interval between time points. e.g. 2 if sampled every 2hrs
  #   OMEGA: 2 * pi / PERIOD. Angular frequency. If omega = 0, matrix Y is 
  #   normalized by time components (not oscillating)
  # 
  # RETURNS:
  #   Y.time.projected: matrix of expression, projected onto the temporal axis.
  
  # track number of genes for dimension purposes
  N.GENES <- nrow(Y)
  # get times vector
  times.vec <- seq(length.out = N.TIMEPTS, by = INTERVAL)
  
  # init output matrix
  Y.time.projected <- matrix(NA, nrow=N.GENES, ncol=N.TISSUES)
  
  # identify row and colnames
  rownames(Y.time.projected) <- rownames(Y)  # same row names
  colnames(Y.time.projected) <- col.names  # from args  
  
  # BEGIN: project onto temporal axis
  for (i in 1:N.TISSUES){
    # get tissue.i across time
    index.start <- (i - 1) * N.TIMEPTS + 1  # 1, 25, 49...
    index.end <- i * N.TIMEPTS
    Y.tissue.i <- Y[, index.start:index.end]  # all genes
    Y.fft <- CalculatePeriodogram(Y.tissue.i, is.matrix=TRUE)
    
    # WHICH FREQUENCY TO EXTRACT? EXTRACT FREQUENCY AT OMEGA.
    freq <- OMEGA / (2 * pi)
    
    # CalculatePeriodogram frequencies are in steps from
    # (0:n.steps/2) by steps 1 / n.steps 
    # we need to divide frequency by INTERVAL to account for sampling interval
    freq.adj <- Y.fft$freq / INTERVAL
    freq.i <- which(freq.adj == freq)
    Y.time.projected[, i] <- as.matrix(Y.fft$p.scaled)[, freq.i]
  }
  return(Y.time.projected)
}

# Test function works -----------------------------------------------------

# test ProjectToPeriodicTIme
#
SHOWTHIS <- FALSE
if (SHOWTHIS == TRUE){
  ROWS <- 12
  COLS <- 288
  n.timepts <- 24  # 24 time points per tissue
  n.tiss <- COLS / n.timepts
  P <- 24  # 24 hour period
  w <- 2 * pi / P
  # w <- 0
  # use a cosine or sine function with period of 24 hours
  t <- seq(from=0, by=2, length.out=COLS)  # 0 to 48 hours sampled every 2 hrs. 12 tissues
  y <- 3 * sin(w * t + pi / 2)
  # y <- rep(1:n.tiss, each=n.timepts)
  Y <- matrix(y, nrow=ROWS, ncol=COLS, byrow = TRUE)
  out.colnames <- make.unique(rep('COL', n.tiss))
  rownames(Y) <- make.unique(rep('ROW', ROWS))
  (Y.t <- ProjectToPeriodicTime(Y, N.TISSUES=n.tiss, N.TIMEPTS=n.timepts, INTERVAL=2, OMEGA=w, out.colnames))
  (Mod(Y.t)) 
  (Y.t)
  
  # DO SVD
  Y.t[6:12] <- 2 * Y.t[6:12]
  s <- svd(Y.t)
  
}

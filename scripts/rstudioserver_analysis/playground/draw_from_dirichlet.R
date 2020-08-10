# Jake Yeung
# Date of Creation: 2020-06-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/playground/draw_from_dirichlet.R
# 

#' ## Thinking of number of unique UMIs and reads per UMI as draws from a Dirichlet distribution

#' generate proportions of different transcripts. We think proportions are heavy tailed. If these proportions come from 
#' a result of a multiplicative product of many independent random variables (maybe that's one way to think of PCR reactions?) 
#' then let's go with the lognormal. But I don't think we need to be very precise about this here. 

Nreads <- 10^8  # total reads sequenced
proportions.pcr <- rlnormrln(n = Nreads, mean = 0, sd = 1)
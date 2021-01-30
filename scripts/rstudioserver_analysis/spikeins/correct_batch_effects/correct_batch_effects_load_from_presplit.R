# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/1-correct_batch_effects_from_peaks_from_LDA.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)
library(parallel)
library(hash)

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()

parser$add_argument('infile', metavar='INFILE',
                    help='Input .rds file contaiing impute_long')
parser$add_argument('outfile', metavar='OUTFILE',
                    help='Output mat adjusted')
parser$add_argument('-ncores', metavar='OUTFILE', type='integer', default = 4,
                    help='Number of cores, default 4')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) {
  print("Arguments:")
  print(args)
}


# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# 
# ncores <- 4
ncores <- args$ncores
# hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load presplit -----------------------------------------------------------

imputed.long <- readRDS(args$infile)

rnames.all <- as.character(unique(imputed.long$rname))
names(rnames.all) <- rnames.all


# Set up environment ------------------------------------------------------

myenv <- new.env()

jsplit <- split(x = imputed.long, f = imputed.long$rname)

print("Loading split data into environment...")
for (jrname in rnames.all){
  myenv[[jrname]] <- jsplit[[jrname]]
}

print("Loading split data into environment...done")

print("Testing one")
print(head(myenv[[jrname]]))

# Clean up  ---------------------------------------------------------------

print("Cleaning up before mclapply")
rm(jsplit)
rm(imputed.long)
gc()

# Correct batch  ----------------------------------------------------------

print(paste("Correcting batch multicore", ncores))
system.time(
  dat.adj.lst <- mclapply(rnames.all, function(jrname){
    jdat <- myenv[[jrname]]
    dat.adj <- jdat %>%
      group_by(rname) %>%
      do(AdjustBatchEffect(.))
    return(dat.adj)
  }, mc.cores = ncores)
)

print(paste("Correcting batch multicore... done"))

dat.adj.long <- dat.adj.lst %>%
  bind_rows()

# dat.adj.lst2 <- lapply(dat.adj.lst, function(jdat){
#   subset(jdat, select = c(rname, cell, log2exprs, cluster, batch, jrep, jrep2, plateadj2, clstradj2, log2exprsadj))
# })

mat.adj <- dat.adj.long %>%
  data.table::dcast(., formula = rname ~ cell, value.var = "log2exprsadj")
# mat.adj.lst <- lapply(dat.adj.lst2, function(dat.adj){
#   mat.adj <- data.table::dcast(dat.adj, formula = rname ~ cell, value.var = "log2exprsadj")
# })

save(mat.adj, file = args$outfile)





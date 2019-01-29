# Jake Yeung
# mat_to_long_mat.R
# Make mat to long mat for covbed 
# 2016-08-03

start <- Sys.time()

library(dplyr)
library(methods)
library(reshape2)

sessionInfo()

# Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

inf <- args[[1]]  # /scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/closestbed_multiple_genes.genomewide.merged/closestbed_multiple_genes.genomewide.merged.bed
outf <- args[[2]]  # /scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/closestbed_multiple_genes.genomewide.merged/closestbed_multiple_genes.genomewide.merged.sqlite3
cnames.str <- args[[3]]  # comma separated
cnames <- strsplit(cnames.str, ",")[[1]]
print("Colnames assigned to table")
print(cnames)
classes.str <- args[[4]]
classes <- strsplit(classes.str, ",")[[1]]

print("Classes assigned to columns")
print(classes)

if (!file.exists(inf)) stop(paste("File does not exist: ", inf))

# Load -------------------------------------------------------------------

mat <- read.table(inf, header = FALSE, sep = "\t", 
  colClasses = classes, 
  col.names = cnames) %>%
  as.data.frame() %>% 
  tbl_df()
print(cnames)

print("Done reading big matrix")
print(Sys.time() - start)

# melt mat
mat.long <- melt(mat, id.vars = c("chromo", "start", "end", "gene", "dist"), variable.name = "sample", value.name = "counts") %>%
  tbl_df()
print("Done melting")
print(head(mat.long))

# print to output
write.table(mat.long, file = outf, append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

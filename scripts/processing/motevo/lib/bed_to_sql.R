# Jake Yeung
# bed_to_mysql.R
#  
# 2016-07-30
# 2019-02-04 Redo

start <- Sys.time()

library(dplyr)
library(methods)
library(reshape2)

sessionInfo()

# Functions ---------------------------------------------------------------

str.to.vector <- function(s, jsep = ","){
  # comma separated string to vector
  # e.g., ZT2,ZT4 -> c("ZT2", "ZT4")
  s.vec <- strsplit(s, jsep)[[1]]
  if (length(s.vec) == 0){
    print(s.vec)
    warning("Length of colnames must be > 0")
  } 
  return(s.vec)
}

str.to.list <- function(s, vec.sep = "-", list.sep = ","){
  # e.g., "chromo-start-end,time" -> list(c(chromo, start, end), time)
  lst <- strsplit(s, list.sep)[[1]]
  lst <- lapply(lst, str.to.vector, vec.sep)
  return(lst)
}

handle.na <- function(s){
  if (s == "NA"){
    return(NA)
  } else {
    return(s)
  }
}


# Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

inf <- args[[1]]  # /scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/closestbed_multiple_genes.genomewide.merged/closestbed_multiple_genes.genomewide.merged.bed
outf <- args[[2]]  # /scratch/el/monthly/jyeung/motevo_dhs_outputs/motevo_outputs/closestbed_multiple_genes.genomewide.merged/closestbed_multiple_genes.genomewide.merged.sqlite3
cnames.str <- handle.na(args[[3]])  # comma separated
if (!is.na(cnames.str)){
  cnames <- str.to.vector(cnames.str)
  has.header <- FALSE
} else {
  has.header <- TRUE
  cnames <- NULL
}

classes.str <- handle.na(args[[4]])
if (!is.na(classes.str)){
  classes <- str.to.vector(classes.str)
}

tbl.name <- args[[5]]
# no dots in tbl.name
tbl.name <- gsub("\\.", "_", tbl.name)

indx.str <- args[[6]]
indx <- str.to.list(indx.str)

if (!file.exists(inf)) stop(paste("File does not exist: ", inf))

# Print arguments
print(paste("Inf:", inf))
print(paste("Outf:", outf))
print(paste("Cnames:")); print(cnames.str)
print(paste("Classes:")); print(classes.str)
print(paste("Table name:", tbl.name))
print(paste("Indx:")); print(indx)

# Main -------------------------------------------------------------------

# my_db <- src_sqlite(outf, create = T)
my_db <- DBI::dbConnect(RSQLite::SQLite(), dbname = outf)

# tab5rows <- read.table(inf, header = has.header, sep = "\t", nrows = 5,
#   colClasses = classes,
#   col.names = cnames) %>%
#   as.data.frame() %>% 
#   tbl_df()
# if (is.na(classes.str)){
#   tab5rows <- read.table(inf, header = has.header, sep = "\t", nrows = 5,
#     col.names = cnames) %>%
#     as.data.frame() %>% 
#     tbl_df()
#   classes <- sapply(tab5rows, class)
#   print("Auto assigning classes:")
#   print(classes)
# }
# 
# # test on small dataset make sure it works
# print("Testing on 5 rows...")
# test.tbl <- copy_to(my_db, tab5rows, "test5rows", temporary = TRUE, indexes = indx)
# print("Done testing 5 rows")
# print(src_tbls(my_db))

# gc()

# load matrix now (could be big!)
if (!has.header){
	mat <- data.table::fread(inf, header = has.header, sep = "\t", 
						colClasses = classes, 
						col.names = cnames) %>%
	  as.data.frame() %>% 
	  tbl_df()
} else {
	mat <- data.table::fread(inf, header = has.header, sep = "\t", 
						colClasses = classes)%>%
	  as.data.frame() %>% 
	  tbl_df()
}

print("Done reading big matrix")
print(Sys.time() - start)

# https://blog.rstudio.com/2017/06/27/dbplyr-1-1-0/
# DBI::dbWriteTable(my_db, "H3K4me1_motevo", mat)

mat.tbl <- copy_to(my_db, mat, tbl.name, temporary = FALSE, indexes = indx)

print("Done copying database")
print(Sys.time() - start)


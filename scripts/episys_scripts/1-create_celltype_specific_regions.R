#!/bin/sh
# Jake Yeung
# create_celltype_specific_regions.sh
#  
# 2019-06-25

library(JFuncs)
library(data.table)

region1="chr7:114,972,165-116,898,716"
region2="chr7:103,325,741-104,425,479"
region3="chr3:90,523,036-90,870,443"
region4="chr11:44,114,099-45,269,522"

regions.lst <- list(Sox6=region1, Hbb=region2, S100a=region3, Ebf1=region4)
regions.lst <- lapply(regions.lst, function(r){
  return(gsub(",", "", r))
})

chromos <- unlist(lapply(regions.lst, GetChromo))
starts <- unlist(lapply(regions.lst, GetStart))
ends <- unlist(lapply(regions.lst, GetEnd))

dat.long <- data.frame(chromo = chromos, start = starts, end = ends)

print(dat.long)
fwrite(dat.long, file = "/Users/yeung/data/scchic/for_episys/regions/regions_to_filter.txt", col.names = FALSE, sep = "\t")

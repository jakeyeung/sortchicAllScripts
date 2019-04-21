# Jake Yeung
# Date of Creation: 2019-04-15
# File: ~/projects/scchic/scripts/Rfunctions/AnnotationFunctions.R
# Annotation functions


AnnotateRegionToGene <- function(peaks, chr.exclude = c("chr20", "chr21")){
  # annotate regions?
  require(JFuncs)
  require(ChIPseeker)
  require(GenomicRanges)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  require(org.Mm.eg.db)
  regions <- data.frame(seqnames = sapply(peaks, GetChromo),
                        start = sapply(peaks, GetStart),
                        end = sapply(peaks, GetEnd),
                        stringsAsFactors = FALSE)
  assertthat::assert_that(nrow(regions) > 0)
  rownames(regions) <- peaks
  regions <- subset(regions, !seqnames %in% chr.exclude)
  
  regions.range <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(ChIPseeker::annotatePeak(regions.range,
                                                  TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                                                  annoDb='org.Mm.eg.db'))
  regions.annotated$region_coord <- names(regions.range)
  return(regions.annotated)
}


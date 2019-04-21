# Jake Yeung
# 1-annotate_bed_to_gene_and_distance.R
# Annotate bed to gene distance
# 2019-01-29

suppressPackageStartupMessages(library("argparse"))
# create parser object
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(annotate))

source("scripts/Rfunctions/Aux.R")

EntrezToSymbol.csv <- function(x, jsplit = ";"){
    # csv of entrez to symbol
    # 59014;76187;76187;72481;72481;76982
    xvec <- strsplit(x, split=jsplit)[[1]]
    xvec.symbol <- getSYMBOL(xvec, data = "org.Mm.eg.db")
    xsymbol <- paste(xvec.symbol, collapse = "@")
    return(xsymbol)
}

parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('infile', metavar='INFILE',
                                            help='Infile bed from hiddenDomains')
parser$add_argument('outfile', metavar='OUTFILE',
                                            help='Outfile with gene name and distance in columns 4 and 5')
parser$add_argument("-m", "--multi", action="store_true", default=FALSE,
                        help="Assign to multiple genes, requires --dist option")
parser$add_argument("-d", "--dist", metavar='INT', default=50000, type="integer",
                        help="Maximum distance from gene for assigning")
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

# load bed file 
bed <- data.table::fread(args$infile)
if (ncol(bed) == 4){
  colnames(bed) = c("chromo", "start", "end", "peakname")
} else {
  colnames(bed) = c("chromo", "start", "end")
}
print(head(bed))
# regions <- data.frame(seqnames = sapply(colnames(bed), GetChromo),
#                       start = sapply(colnames(bed), GetStart),
#                       end = sapply(colnames(bed), GetEnd), 
#                       stringsAsFactors = FALSE)
regions <- data.frame(seqnames = bed$chromo,
                      start = bed$start,
                      end = bed$end,
                      stringsAsFactors = FALSE)
rownames(regions) <- regions$peakname
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))  # remove X and Y?

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))

if (!args$multi){
    regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                    TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                    annoDb='org.Mm.eg.db'))
} else {
    regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                    TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                    addFlankGeneInfo=TRUE,
                                                    flankDistance=args$dist,
                                                    annoDb='org.Mm.eg.db'))

}
regions.annotated$region_coord <- names(regions.range)

print(head(regions.annotated))
print(colnames(regions.annotated))

# keep gene name and distance
if (!args$multi){
    regions.out <- subset(regions.annotated, select = c(seqnames, start, end, SYMBOL, distanceToTSS))
} else {
    # keep column names the same
    regions.out <- subset(regions.annotated, select = c(seqnames, start, end, flank_geneIds, flank_gene_distances))
    # replace flank_geneIDs from entrezID to symbol
    regions.out$flank_geneIds <- sapply(regions.out$flank_geneIds, function(x){
                                        if (!is.na(x)){
                                          xsym <- EntrezToSymbol.csv(x, jsplit = ";")
                                        } else {
                                          xsym <- x
                                        }
                      })
    # use @ as separation to fit with /home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/convert_compressed_bed_to_long.py
    regions.out$flank_gene_distances <- sapply(regions.out$flank_gene_distances, function(x) gsub(";", "@", x))
}

print(regions.out)


# write to table
data.table::fwrite(regions.out, file = args$outfile, append=FALSE, quote=FALSE, sep = "\t", col.names=FALSE, na = "NA")

# Jake Yeung
# 3-pwm_to_seqlogo.R
# https://davetang.org/muse/2013/01/30/sequence-logos-with-r/
# 2017-06-14

library(seqLogo)

proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

args <- commandArgs(trailingOnly=TRUE)
# inf <- "/scratch/el/monthly/jyeung/databases/WMs/SwissRegulon_mm10/split_by_motif/Arntl.pwm"
inf <- args[[1]]
outf <- args[[2]]
if (length(args) == 3){
  outformat <- args[[3]]  # pdf or jpeg
} else {
  outformat <- "jpg"  # default saves space 
}

print(paste("outformat:", outformat))

ncol <- max(count.fields(inf, sep = ""))
data <- read.table(inf, sep="", fill=TRUE, col.names=1:ncol, stringsAsFactors=FALSE)

motif <- data[2, 2]
print(motif)

colnames(data) <- tolower(data[3, 1:ncol(data)])
data <- data[4:(nrow(data) - 1), 2:5]
# print(data)
data <- as.matrix(data, ncols = 4)
data <- matrix(as.numeric(data), ncol = 4)
# print(data)

pwm <- apply(data, 1, proportion)
pwm <- makePWM(pwm)

if (outformat == "pdf"){
  pdf(outf)

} else if (outformat == "jpeg" | outformat == "jpg"){
  jpeg(outf, quality = 30)

} else {
  stop(paste("Outformat must be pdf or jpg", outformat))
}
seqLogo(pwm)
dev.off()

print(paste("Saved to:", outf))


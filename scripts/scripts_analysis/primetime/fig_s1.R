# Jake Yeung
# Date of Creation: 2019-02-21
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/fig_s1.R
# Fig S1

setwd("~/Dropbox/scCHiC_figs/figures/FIG_S1")

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K9me3 = "red1", H3K27me3 = "darkorange1")

q1<-read.table('PZ-K562-G1-H3K4me1.umi')
q2<-read.table('PZ-K562-G1-H3K4me3.umi')
q3<-read.table('PZ-K562-G1-H3K9me3.umi')
q4<-read.table('PZ-K562-G1-H3K27me3.umi')
q4<-q4[q4[,4]<1e7,]

h1<-hist(log10(q1[,4]),breaks=100)
h2<-hist(log10(q2[,4]),breaks=100)
h3<-hist(log10(q3[,4]),breaks=100)
h4<-hist(log10(q4[,4]),breaks=100)


pdf('Unique_Cuts.pdf', paper = "a4", useDingbats = FALSE)
# png('Unique_Cuts.png', width=500, height = 1000)

ncells <- paste(sapply(list(q1, q2, q3, q4), function(x) nrow(x)), collapse = ",")
jmain <- paste("K562, N=", ncells)
par(mfrow=c(4,1),mar=c(0,0,0,0))
a=0.8
b=5.2
plot(h1$mids,h1$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[1]],lwd=4,main=jmain)
points(c(a,b),c(0,0),type='l',col=colsvec[[1]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')
plot(h2$mids,h2$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[2]],lwd=4)
points(c(a,b),c(0,0),type='l',col=colsvec[[2]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')
plot(h3$mids,h3$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[3]],lwd=4)
points(c(a,b),c(0,0),type='l',col=colsvec[[3]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')
plot(h4$mids,h4$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[4]],lwd=4)
axis(1, at=1:5, labels=c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)))
points(c(a,b),c(0,0),type='l',col=colsvec[[4]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')



q1<-read.table('PZ-BM-H3K4me1.umi')
q2<-read.table('PZ-BM-H3K4me3.umi')
q3<-read.table('PZ-BM-H3K9me3.umi')
q4<-read.table('PZ-BM-H3K27me3.umi')
q4<-q4[q4[,4]<1e7,]

# h1<-hist(log10(q1[,4]),breaks=100)
# h2<-hist(log10(q2[,4]),breaks=100)
# h3<-hist(log10(q3[,4]),breaks=100)
# h4<-hist(log10(q4[,4]),breaks=100)

par(mfrow=c(4,1),mar=c(0,0,0,1))
a=0.8
b=5.2

ncells <- paste(sapply(list(q1, q2, q3, q4), function(x) nrow(x)), collapse = ",")
jmain <- paste("BM, N=", ncells)

plot(h1$mids,h1$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[1]],lwd=4, main = jmain)
points(c(a,b),c(0,0),type='l',col=colsvec[[1]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')
plot(h2$mids,h2$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[2]],lwd=4)
points(c(a,b),c(0,0),type='l',col=colsvec[[2]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')
plot(h3$mids,h3$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[3]],lwd=4)
points(c(a,b),c(0,0),type='l',col=colsvec[[3]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')
plot(h4$mids,h4$counts,xlim=c(a,b),type='h',axes=F,col=colsvec[[4]],lwd=4)
axis(1, at=1:5, labels=c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)))
points(c(a,b),c(0,0),type='l',col=colsvec[[4]],lwd=4)
grid(ny=0,lty=1,col='darkgrey')



dev.off()

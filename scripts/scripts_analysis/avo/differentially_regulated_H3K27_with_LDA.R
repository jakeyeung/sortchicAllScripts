# rm(list=ls())
# setwd("~/Dropbox/scCHiC_figs/FIG4_BM/")

xmin<-95.3e6
xmax<-97.3e6
ymin=-0.1
ymax=.5
chr<-'c11_'

load("~/Dropbox/scCHiC_figs/FIG4_BM/H3K27me3.datadir_mc_f.Rda")
Q<-object@mc_fp
mc_index<-object@mc
mc_colors<-object@colors
Q<-data.frame(Q)
rownames(Q)<-paste('c',rownames(Q),sep='')
for (i in 1:ncol(Q)){
	colnames(Q)[i]<-paste('H3K4me27_','mc',i,sep='')
}

allgroup <- seq(11)
ingroup <- c(5,7)  # 5 is most outer group, 7 is less    # group 5 shows 0.6 fold change. group 7 shows 0.3 FC
outgroup <- allgroup[!allgroup %in% ingroup]

if (length(ingroup) > 1){
  Q<-data.frame(cbind(rowMeans(Q[,c(ingroup)]),rowMeans(Q[,c(outgroup)])))
} else {
  Q<-data.frame(cbind(Q[, ingroup],rowMeans(Q[,c(outgroup)])))
}

q<-strsplit(rownames(Q),'_')
x<-NULL
for (i in 1:nrow(Q)){
	x[i]<-as.integer(q[[i]][2])
}
Q$coord<- x+50000
x<-Q$coord[grep(chr,rownames(Q))]

plot(log2(Q[,1]),log2(Q[,2]),pch=16,cex=.25,xlab='average log2-fold change metacell 5 & 7',ylab='average log2-fold change other metacells')

#HoxA chr6 midpoint 52.2 Mb
x0<-52200000
left<-x0-100000
right<-x0+100000
Qchr<-Q[grep('c6_',rownames(Q)),]
Qchr<-Qchr[(Qchr$coord>left & Qchr$coord<right),]
Qchr<-Qchr[!is.na(Qchr$coord),]
points(log2(Qchr[,1]),log2(Qchr[,2]),pch=16,cex=1,col='pink')

#HoxB chr11 midpoint 96.3 Mb
x0<-96300000
left<-x0-100000
right<-x0+100000
Qchr<-Q[grep('c11_',rownames(Q)),]
Qchr<-Qchr[(Qchr$coord>left & Qchr$coord<right),]
Qchr<-Qchr[!is.na(Qchr$coord),]
points(log2(Qchr[,1]),log2(Qchr[,2]),pch=16,cex=1,col='blue')

#HoxC chr15 midpoint 103 Mb
x0<-103000000
left<-x0-100000
right<-x0+100000
Qchr<-Q[grep('c15_',rownames(Q)),]
Qchr<-Qchr[(Qchr$coord>left & Qchr$coord<right),]
Qchr<-Qchr[!is.na(Qchr$coord),]
points(log2(Qchr[,1]),log2(Qchr[,2]),pch=16,cex=1,col='green')

#HoxC chr2 midpoint 74.7 Mb
x0<-74700000
left<-x0-100000
right<-x0+100000
Qchr<-Q[grep('c2_',rownames(Q)),]
Qchr<-Qchr[(Qchr$coord>left & Qchr$coord<right),]
Qchr<-Qchr[!is.na(Qchr$coord),]
points(log2(Qchr[,1]),log2(Qchr[,2]),pch=16,cex=1,col='cyan')
legend(0.36,0.036, c('HoxA','HoxB','HoxC','HoxD'),col=c('red','blue','green','cyan'),pch=c(16,16,16,16),bty='n')

# where are metacells 5 and 7 on the KNN graph?

# top right two clusters

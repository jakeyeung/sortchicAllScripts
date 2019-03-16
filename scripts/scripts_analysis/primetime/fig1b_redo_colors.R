# Jake Yeung
# Date of Creation: 2019-02-21
# File: ~/projects/scchic/scripts/scripts_analysis/primetime/fig1b_redo_colors.R
# Redo colors for Figure S1B


setwd("~/Dropbox/scCHiC_figs/figures/FIG1")


# Set constants -----------------------------------------------------------


jlims <- list(c(15e6, 115e6), c(35e6, 80e6), c(50e6, 60e6))


# Load and plot -----------------------------------------------------------



qq1<-read.table('K562-G1-H3K4me1-100kb.txt')
q1<-qq1[rowSums(qq1)<4000,]
q1<-q1[,colSums(q1)>10000]

qq2<-read.table('K562-G1-H3K27me3-100kb.txt')
q2<-qq2[rowSums(qq2)<4000,]
q2<-q2[,colSums(q2)>10000 & colSums(q2)<1e7,]

qq3<-read.table('K562-G1-H3K4me3-100kb.txt')
q3<-qq3[rowSums(qq3)<4000,]
q3<-q3[,colSums(q3)>10000]

qq4<-read.table('K562-G1-H3K9me3-100kb.txt')
q4<-qq4[rowSums(qq4)<4000,]
q4<-q4[,colSums(q4)>10000 & colSums(q4)<1e7,]

qq5<-read.table('ChIP_data_binned.dat')
q5<-qq5[qq5$V1=="chr1",]
q5$coord<-q5$V2+50000

q1<-q1[,order(-colSums(q1))]
q2<-q2[,order(-colSums(q2))]
q3<-q3[,order(-colSums(q3))]
q4<-q4[,order(-colSums(q4))]

rownames(q1)<-paste('c',rownames(q1),sep='')
rownames(q2)<-paste('c',rownames(q2),sep='')
rownames(q3)<-paste('c',rownames(q3),sep='')
rownames(q4)<-paste('c',rownames(q4),sep='')

q1<-q1[grep('c1_',rownames(q1)),]
q2<-q2[grep('c1_',rownames(q2)),]
q3<-q3[grep('c1_',rownames(q3)),]
q4<-q4[grep('c1_',rownames(q4)),]

Q1<-strsplit(rownames(q1),'_')
x1<-NULL
for (i in 1:nrow(q1)){
  x1[i]<-as.integer(Q1[[i]][2])
}
q1$coord<-x1+50000

Q2<-strsplit(rownames(q2),'_')
x2<-NULL
for (i in 1:nrow(q2)){
  x2[i]<-as.integer(Q2[[i]][2])
}
q2$coord<-x2+50000

Q3<-strsplit(rownames(q3),'_')
x3<-NULL
for (i in 1:nrow(q3)){
  x3[i]<-as.integer(Q3[[i]][2])
}
q3$coord<-x3+50000

Q4<-strsplit(rownames(q4),'_')
x4<-NULL
for (i in 1:nrow(q4)){
  x4[i]<-as.integer(Q4[[i]][2])
}
q4$coord<-x4+50000

SD<-apply(q4,2,sd)
M<-apply(q4,2,mean)
SD<-SD[1:(length(M)-1)]
M<-M[1:(length(M)-1)]

# png(width=600,height=1000)

colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K9me3 = "red1", H3K27me3 = "darkorange1")
pdf("~/Dropbox/scCHiC_figs/figures/FIG1/Fig1_peaks_recolored.pdf", useDingbats = FALSE)
# png("~/Dropbox/scCHiC_figs/figures/FIG1/Fig1_peaks_recolored.pdf", width = 600, height = 1000, useDingbats = FALSE)

for (jlim in jlims){
  # a=15e6
  # b=115e6
  a <- as.integer(jlim[[1]])
  b <- as.integer(jlim[[2]])
  N=3
  L=5
  select1<-c(4,6,9)
  select2<-c(18,32,42)
  select3<-c(4,7,56)
  select4<-c(14,37,38)
  
  
  par(mfrow=c(4*N+4,1),mar=c(0,0,1,0))
  
  plot(q5$coord,q5$H3K4me1_ChIP/q5$input_ChIP,pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K4me1,ylim=c(0,6),lwd=L, main = paste0(a,"-", b))
  for (i in 1:N){plot(q1$coord,q1[,select1[i]],pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K4me1,lwd=L)}
  
  plot(q5$coord,q5$H3K4me3_ChIP/q5$input_ChIP,pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K4me3,ylim=c(0,6),lwd=L)
  for (i in 1:N){plot(q3$coord,q3[,select3[i]],pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K4me3,lwd=L)}
  
  y<-q5$H3K9me3_ChIP/q5$input_ChIP
  y[y>4.5]<-0
  plot(q5$coord,y,pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K9me3,ylim=c(1.25,5),lwd=L)
  for (i in 1:N){plot(q4$coord,q4[,select4[i]],pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K9me3,lwd=L)}
  
  plot(q5$coord,q5$H3K27me3_ChIP/q5$input_ChIP,pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K27me3,ylim=c(0,5),lwd=L)
  for (i in 1:N){plot(q2$coord,q2[,select2[i]],pch=16,cex=.5,xlim=c(a,b),type='h',axes=F,col=colsvec$H3K27me3,lwd=L)}
  
}

dev.off()


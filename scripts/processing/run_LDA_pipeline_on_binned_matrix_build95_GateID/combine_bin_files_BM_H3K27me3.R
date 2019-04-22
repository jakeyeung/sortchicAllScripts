file2<-"/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m1-H3K27me3-2_AH3VGVBGX9_S2/cuts.bin"
file1<-"/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14/cuts.bin"
file4<-"/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m2-H3K27me3-2_H2GV2BGX9_S18/cuts.bin"
file3<-"/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed/PZ-BM-m2-H3K27me3-1_AH3VGVBGX9_S6/cuts.bin"

# fileout<-"/hpc/hub_oudenaarden/avo/scChiC/metacell/BM-H3K27me3-buildnew.txt"
fileout<-"/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_tables_raw/BM-H3K27me3-build95.txt"

q1<-read.table(file1)
q2<-read.table(file2)
q3<-read.table(file3)
q4<-read.table(file4)

rownames(q1)<-paste(as.integer(q1[,1]),as.integer(q1[,2]),as.integer(q1[,3]),sep="_")
rownames(q2)<-paste(as.integer(q2[,1]),as.integer(q2[,2]),as.integer(q2[,3]),sep="_")
rownames(q3)<-paste(as.integer(q3[,1]),as.integer(q3[,2]),as.integer(q3[,3]),sep="_")
rownames(q4)<-paste(as.integer(q4[,1]),as.integer(q4[,2]),as.integer(q4[,3]),sep="_")


q1[,1:3]<-NULL
q2[,1:3]<-NULL
q3[,1:3]<-NULL
q4[,1:3]<-NULL


for (i in 1:384){
colnames(q1)[i]<-paste('BM_H3K27me3_m1_rep1_cell',i,sep="")
colnames(q2)[i]<-paste('BM_H3K27me3_m1_rep2_cell',i,sep="")
colnames(q3)[i]<-paste('BM_H3K27me3_m2_rep1_cell',i,sep="")
colnames(q4)[i]<-paste('BM_H3K27me3_m2_rep2_cell',i,sep="")}


q<-data.frame(cbind(q1,q2,q3,q4))

q<-q[rowSums(q)>1,]
write.table(q,fileout,quote=FALSE,sep='\t')

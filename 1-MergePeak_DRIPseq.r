#!/apps/R-4.0.2/bin/R

library(Rsamtools)
library(DESeq2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

# Step1: Process and obtain data

case1="../DRIP_shTDP43_1/DRIP_shTDP43_1_peaks.narrowPeak"
case2="../DRIP_shTDP43_2/DRIP_shTDP43_2_peaks.narrowPeak"
case3="../DRIP_shTDP43_3/DRIP_shTDP43_3_peaks.narrowPeak"
control1="../DRIP_shNC_1/DRIP_shNC_1_peaks.narrowPeak"
control2="../DRIP_shNC_2/DRIP_shNC_2_peaks.narrowPeak"
control3="../DRIP_shNC_3/DRIP_shNC_3_peaks.narrowPeak"

case1.gr=import(case1)
case2.gr=import(case2)
case3.gr=import(case3)
control1.gr=import(control1)
control2.gr=import(control2)
control3.gr=import(control3)

case.peak.gr<-Reduce(union, list(
  case1.gr,
  case2.gr,
  case3.gr))

control.peak.gr=Reduce(union, list(
  control1.gr,
  control2.gr,
  control3.gr))

case.peak<-as.data.frame(case.peak.gr)
control.peak<-as.data.frame(control.peak.gr)
write.table(case.peak,file="case.peak.bed",row.names=F,col.names=F,sep="\t",quote=F)
write.table(control.peak,file="control.peak.bed",row.names=F,col.names=F,sep="\t",quote=F)

all.peak.gr=Reduce(union, list(
  case1.gr,
  case2.gr,
  case3.gr,
  control1.gr,
  control2.gr,
  control3.gr))

all.peak<-as.data.frame(all.peak.gr)
write.table(all.peak,file="all.peak.bed",row.names=F,col.names=F,sep="\t",quote=F)

case1="../sortedbed/DRIP_shTDP43_1.sorted.bed"
case2="../sortedbed/DRIP_shTDP43_2.sorted.bed"
case3="../sortedbed/DRIP_shTDP43_3.sorted.bed"
control1="../sortedbed/DRIP_shNC_1.sorted.bed"
control2="../sortedbed/DRIP_shNC_2.sorted.bed"
control3="../sortedbed/DRIP_shNC_3.sorted.bed"

case1.gr=import(case1)
case2.gr=import(case2)
case3.gr=import(case3)
control1.gr=import(control1)
control2.gr=import(control2)
control3.gr=import(control3)

case1c=countOverlaps(all.peak.gr, case1.gr)
case2c=countOverlaps(all.peak.gr, case2.gr)
case3c=countOverlaps(all.peak.gr, case3.gr)
control1c=countOverlaps(all.peak.gr, control1.gr)
control2c=countOverlaps(all.peak.gr, control2.gr)
control3c=countOverlaps(all.peak.gr, control3.gr)

mat=data.frame(case1c, case2c, case3c,control1c, control2c, control3c)

sf=apply(mat,2,sum)
sf=sf/min(sf)
matnorm=sweep(mat,2,sf,FUN="/")

save(mat,file="mat.rda")

caseread=apply(matnorm[,1:3],1,mean)
controlread=apply(matnorm[,4:6],1,mean)

MAT<-mat

#mat<-mat[which(
#               (matnorm$control1c>20 & matnorm$control2c>20 & matnorm$control3c>20 & matnorm$control4c>20) |
#               (matnorm$case1c>20 & matnorm$case2c>20 & matnorm$case3c>20 & matnorm$case3c>20)),]

mat<-mat[(mat$control1c>20 & mat$control2c>20 & mat$control3c>20) |
               (mat$case1c>20 & mat$case2c>20 & mat$case3c>20),]

condition=c(rep("treated",3),rep("control",3))
colData=data.frame(condition=condition)
dds=DESeqDataSetFromMatrix(countData =mat,colData = colData,design= ~condition)
dds=DESeq(dds)
result=results(dds)
result<-as.data.frame(result)

id=as.numeric(rownames(result))
resultmat=MAT[id,]
caseread=caseread[id]
controlread=controlread[id]
fc=(caseread+0.001)/(controlread+0.001)
all.peak.gr=as.data.frame(all.peak.gr)
all.peak.gr=all.peak.gr[id,1:3]
allmat=data.frame(all.peak.gr=all.peak.gr,controlread=controlread,caseread=caseread,fc=fc,result=result)

allmat<-allmat[order(allmat$result.pval),]
allmat$fdr.self<-p.adjust(allmat$result.pval,method="fdr")
allmatFDR<-allmat[allmat$fdr.self<0.1,]
allmatFDR<-allmatFDR[!is.na(allmatFDR$fdr.self),]

write.table(allmat,file="DESeq.xls",sep="\t",quote=FALSE,row.names=FALSE)
write.table(allmatFDR,file="DESeq.FDR0.1.xls",sep="\t",quote=FALSE,row.names=FALSE)



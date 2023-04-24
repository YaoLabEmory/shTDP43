data<-read.table(file="DownInteraction.Enhancer.Promoter.bed.EnhancerlogFC.PromoterlogFC.genelogFC.InteractionlogFC.Bin1RlooplogFC.Bin2RlooplogFC.5hmCEnhancerlogFC.GeneFDR")
data<-data[,c(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)]
names(data)<-c("EnhancerChr","EnhancerStart","EnhancerEnd","PromoterChr","PromoterStart","PromoterEnd","EnhancerlogFC","PromoterlogFC","Gene","GenelogFC","InteractionlogFC","Bin1logFC","Bin2logFC","hmCEnhancerlogFC","GeneFDR")
data<-data[!is.na(data$GenelogFC),]

#data<-data[data$GeneFDR<0.05,]

DownEnhDown5hmC<-unique(data[data$EnhancerlogFC<0 & data$hmCEnhancerlogFC<0,]$Gene)
write.table(DownEnhDown5hmC,file="DownEnhancerDown5hmCGenelist",row.names=F,col.names=T,sep="\t",quote=F)

DownEnhDownGene<-data[data$EnhancerlogFC<0 & data$GenelogFC<0,]
DownEnhDownGenesort<-DownEnhDownGene[order(DownEnhDownGene$GenelogFC),]
dim(DownEnhDownGenesort)
DownEnhDownGeneSsortMoreChange<-head(DownEnhDownGenesort,189)
DownEnhDownGeneSsortLessChange<-tail(DownEnhDownGenesort,189)

DownEnhDownGeneSsortMoreChange$Type="Big"
DownEnhDownGeneSsortLessChange$Type="Small"

t.test(DownEnhDownGeneSsortMoreChange$EnhancerlogFC,DownEnhDownGeneSsortLessChange$EnhancerlogFC)

df<-rbind(DownEnhDownGeneSsortMoreChange,DownEnhDownGeneSsortLessChange)
write.table(df,file="DownEnhancerDownGeneMoreOrLessChange.xls",row.names=F,col.names=T,sep="\t",quote=F)

DownEnhDown5hmCDownGene<-data[data$EnhancerlogFC<0 & data$GenelogFC<0 & data$hmCEnhancerlogFC<0,]
DownEnhDown5hmCDownGeneSsort<-DownEnhDown5hmCDownGene[order(DownEnhDown5hmCDownGene$GenelogFC),]
dim(DownEnhDown5hmCDownGeneSsort)
DownEnhDown5hmCDownGeneSsortMoreChange<-head(DownEnhDown5hmCDownGeneSsort,126)
DownEnhDown5hmCDownGeneSsortLessChange<-tail(DownEnhDown5hmCDownGeneSsort,126)

DownEnhDown5hmCDownGeneSsortMoreChange$Type="Big"
DownEnhDown5hmCDownGeneSsortLessChange$Type="Small"

t.test(DownEnhDown5hmCDownGeneSsortMoreChange$hmCEnhancerlogFC,DownEnhDown5hmCDownGeneSsortLessChange$hmCEnhancerlogFC)

data<-rbind(DownEnhDown5hmCDownGeneSsortMoreChange,DownEnhDown5hmCDownGeneSsortLessChange)
write.table(data,file="DownEnhancerDown5hmCDownGeneMoreOrLessChange.xls",row.names=F,col.names=T,sep="\t",quote=F)

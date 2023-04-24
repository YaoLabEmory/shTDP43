GO<-c("response to stimulus (GO:0050896)","response to chemical (GO:0042221)","regulation of growth (GO:0040008)","regulation of developmental growth (GO:0048638)","nervous system development (GO:0007399)")

DDD<-read.table(file="DownDownDown.GO.xls",header=T,quote="",fill=FALSE,sep="\t")[,c(1,3,6,8)]
names(DDD)<-c("GO","DDDGeneNumber","DDDFoldEnrich","DDDFDR")
DnotDDD<-read.table(file="DnotDDD.GO.xls",,header=T,quote="",fill=FALSE,sep="\t")[,c(1,3,6,8)]
names(DnotDDD)<-c("GO","DnotDDDGeneNumber","DnotDDDFoldEnrich","DnotDDDFDR")

data<-merge(DDD,DnotDDD,by="GO")

data<-data[(data$DDDFDR<0.05  | data$DnotDDDFDR<0.05),c("GO","DDDFDR","DnotDDDFDR","DDDGeneNumber")]

data$logDnotDDDFDR<-(-log(data$DnotDDDFDR,10))
data$logDDDFDR<-(-log(data$DDDFDR,10))

Selected<-data[data$GO %in% GO,]
NonSelected<-data[!(data$GO %in% GO),]

Selected$color="blue"
NonSelected$color="grey"

data<-rbind(NonSelected,Selected)
write.table(data,file="CombineGO.xls",row.names=F,col.names=T,sep="\t",quote=F)

library(ggplot2)
library(ggallin)
p<-ggplot(data,aes(logDnotDDDFDR,logDDDFDR,color=color))+
   geom_point(aes(size=DDDGeneNumber),data=data)+
   scale_color_manual(values=c("blue","grey"))+
   theme_classic()+
   theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))+
        guides(colour = guide_legend(override.aes = list(size=6)))+
   theme(axis.text.x = element_text( color="black",
                           size=30, angle=0,hjust=0.5),
          axis.text.y = element_text( color="black",
                           size=30, angle=90,hjust=0.5))+
   theme(legend.text=element_text(size=30))+
   scale_x_continuous(trans = pseudolog10_trans,limits = c(0, 10))+
   scale_y_continuous(trans = pseudolog10_trans,limits = c(0, 10))+
   theme(axis.line = element_line(colour = 'black', size = 1.5))+
   xlab("-logFDR DnotDDD")+
   ylab("-logFDR DDD")+
   geom_vline(xintercept = -log(0.05,10),linetype = "dashed",size = 1)+
   geom_hline(yintercept = -log(0.05,10),linetype = "dashed",size = 1)+
   theme(legend.position="none")

pdf(file="Correlation_DDD_DnotDDD.pdf",6,6)
p
dev.off()

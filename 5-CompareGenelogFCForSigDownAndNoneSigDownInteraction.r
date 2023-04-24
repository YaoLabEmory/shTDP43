data<-read.table(file="AllInteraction.Enhancer.Promoter.bed.EnhancerlogFC.PromoterlogFC.genelogFC.InteractionlogFC.Bin1RlooplogFC.Bin2RlooplogFC.5hmCEnhancerlogFC.GeneFDR.InteractionFDR",sep="\t")
names(data)<-c("Bin1Chr","Bin1Start","Bin1End","Bin2Chr","Bin2Start","Bin2End","EnhancerChr","EnhancerStart","EnhancerEnd","PromoterChr","PromoterStart","PromoterEnd","EnhancerRlooplogFC","PromoterRlooplogFC","Gene","GenelogFC","InteractionlogFC","Bin1RlooplogFC","Bin2RlooplogFC","5hmCEnhancerlogFC","GeneFDR","InteractionFDR")
data<-data[,c(7:9,16,17,21,22)]
data$enhancer<-paste(data$EnhancerChr,data$EnhancerStart,data$EnhancerEnd)
data<-data[,-c(1:3)]
data<-data[!is.na(data$GeneFDR),]
data<-data[!is.na(data$GenelogFC),]
data<-data[is.finite(data$GenelogFC),]
data$type<-"nochange"
data[data$GenelogFC>0,]$type<-"nonsigup"
data[data$GeneFDR<0.05 & data$GenelogFC>0,]$type<-"sigup"
data[data$GenelogFC<0,]$type<-"nonsigdown"
data[data$GeneFDR<0.05 & data$GenelogFC<0,]$type<-"sigdown"

sigdown<-data[data$InteractionFDR<0.05 & data$InteractionlogFC<0,]
other<-data[data$InteractionFDR>0.05 & data$InteractionlogFC<0,]

other<-other[!(other$enhancer %in% sigdown$enhancer),]

dim(sigdown)
dim(other)

table(sigdown$type)
table(other$type)

summary(sigdown$GenelogFC)
summary(other$GenelogFC)

sigdown$InterType<-"sigdown"
other$InterType<-"other"

write.table(sigdown,file="sigdowninteraction.xls",row.names=F,col.names=T,quote=F,sep="\t")
write.table(other,file="otherinteraction.xls",row.names=F,col.names=T,quote=F,sep="\t")

sigdown<-sigdown[sigdown$GenelogFC<0,]
other<-other[other$GenelogFC<0,]

df<-rbind(sigdown,other)

t.test(sigdown$GenelogFC,other$GenelogFC)

library(ggplot2)
p<-ggplot(df, aes(x=InterType, y=GenelogFC, fill=InterType))+
        geom_boxplot(outlier.shape = NA)+
        scale_fill_manual(values = c("sigdown" = "red", "other" = "grey"))+
        scale_colour_manual(values = c("sigdown" = "red", "other" = "grey"))
p<-p+theme()+
            xlab("Interaction Type")+
            ylab("Target Gene logFC")+
            theme_classic()+
            theme( axis.title.x = element_text(size = 0),
            axis.text.y = element_text(size = 30,colour="black"),
            axis.text.x = element_text(size = 15,colour="black"),
            axis.title.y = element_text(size = 30,colour="black"))+
            theme(axis.text.x = element_text(angle = 0,hjust = 0.5))+
            theme(legend.text=element_text(size=30,colour="black"))+
            theme(legend.title=element_text(size=0))+
            theme(legend.position='none')+
            theme(axis.line = element_line(colour = 'black', size = 2))+
            ylim(-0.5,0)

pdf(file="logFC.pdf",6,6)
p
dev.off()

data<-read.table(file="DownInteraction.Enhancer.Promoter.bed.EnhancerlogFC.PromoterlogFC.genelogFC.InteractionlogFC.Bin1RlooplogFC.Bin2RlooplogFC.5hmCEnhancerlogFC.GeneFDR",sep="\t")
names(data)<-c("Bin1Chr","Bin1Start","Bin1End","Bin2Chr","Bin2Start","Bin2End","EnhancerChr","EnhancerStart","EnhancerEnd","PromoterChr","PromoterStart","PromoterEnd","EnhancerRlooplogFC","PromoterRlooplogFC","Gene","GenelogFC","InteractionlogFC","Bin1RlooplogFC","Bin2RlooplogFC","5hmCEnhancerlogFC","GeneFDR")
data<-data[,c(7:9,13,16,17,21)]
data$enhancer<-paste(data$EnhancerChr,data$EnhancerStart,data$EnhancerEnd)
data<-data[,-c(1:3)]
data<-data[!is.na(data$GeneFDR),]
data<-data[!is.na(data$GenelogFC),]
data<-data[is.finite(data$GenelogFC),]
#data<-data[data$GeneFDR<0.05,]
data<-data[data$GenelogFC<0,]
data<-data[data$EnhancerRlooplogFC<0,]

data<-data[order(data$EnhancerRlooplogFC),]

length<-round(length(data$EnhancerRlooplogFC)/10)

length

firdownRloop<-data[1:length,]
#secdownRloop<-data[(length+1):(2*length),]
#thidownRloop<-data[(2*length+1):(3*length),]
lastdownRloop<-data[(9*length+1):(10*length),]

t.test(firdownRloop$GenelogFC,lastdownRloop$GenelogFC)

firdownRloop$Type="firdownRloop"
#secdownRloop$Type="secdownRloop"
#thidownRloop$Type="thidownRloop"
lastdownRloop$Type="lastdownRloop"

df<-rbind(firdownRloop,lastdownRloop)

df$Type<-factor(df$Type,level=c("lastdownRloop","firdownRloop"))

library(ggplot2)
p<-ggplot(df, aes(x=Type, y=GenelogFC, fill=Type))+
        geom_boxplot(outlier.shape = NA)+
        scale_fill_manual(values = c("firdownRloop" = "red","lastdownRloop" = "grey"))+
        scale_colour_manual(values = c("firdownRloop" = "red","lastdownRloop" = "grey"))
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

pdf(file="logFC.Rloop.pdf",6,6)
p
dev.off()

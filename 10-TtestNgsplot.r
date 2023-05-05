library("ggplot2")
library(reshape2)

data<-read.table(file="avgprof.txt",header=T)

data<-data[c(22:89),]

p<-t.test(data[,1],data[,2])$p.value
write.table(p,file="pvalue",row.names=F,col.names=F,sep="\t",quote=F)

data<-melt(data)
names(data)<-c("group","value")

data$group<-factor(data$group,levels=c("shTDP43_DRIP","shNC_DRIP"))



p<-ggplot(data, aes(x=group, y=value, fill=group)) +geom_boxplot(position = position_dodge(),outlier.size = -1)
        #scale_fill_manual(values = c("IDHmut" = "red", "IDHwt" = "blue"))+
        #scale_colour_manual(values = c("IDHmut" = "red", "IDHwt" = "blue"))
p<-p+theme()+
            xlab("")+
            ylab("Normalized Reads")+
            theme_classic()+
            theme( axis.title.x = element_text(size = 0),
            axis.text.y = element_text(size = 20,colour="black"),
            axis.text.x = element_text(size = 20,colour="black"),
            axis.title.y = element_text(size = 20,colour="black"))+
            theme(axis.text.x = element_text(angle = 0,hjust = 0.5))+
            theme(legend.text=element_text(size=20,colour="black"))+
            theme(legend.title=element_text(size=0))+
            theme(legend.position='none')+
#           ylim(0,100)+
            theme(axis.line = element_line(colour = 'black', size = 0.8))+
	    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file="ngspot-boxplot.pdf",3,5)
p
dev.off()

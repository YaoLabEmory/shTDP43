library("GOplot")
david<-read.table(file="ECdavid",header=T,sep="\t",quote="")
genelist<-read.table(file="ECgenelist",header=T,sep="\t")
genes<-read.table(file="ECgenes",header=F)
names(genes)<-"ID"
genes<-merge(genes,genelist,by="ID")[,c(1,2)]
process<-read.table(file="ECprocess",header=F,sep="\t")[,1]

circ <- circle_dat(david, genelist)
chord <- chord_dat(data = circ, genes = genes, process = process)
pdf(file="GOchord.pdf",16.5,18.2)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 10,process.label=10)
dev.off()


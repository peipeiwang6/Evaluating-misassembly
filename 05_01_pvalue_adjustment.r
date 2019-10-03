args=commandArgs(TRUE)
infile = args[1]
outfile = args[2]
infile <- read.table(infile,stringsAsFactors=F,sep="\t",header=F)
pval<-p.adjust(infile[,5],'fdr')
infile<-cbind(infile,pval)
write.table(infile, file = outfile,row.names=F,col.names=F,quote=F,sep="\t")

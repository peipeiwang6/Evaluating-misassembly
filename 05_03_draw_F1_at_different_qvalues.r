args=commandArgs(TRUE)
infile = args[1] ### F-measure between results of two CNVnator runs

library(ggplot2)
library(grid)
setwd("D:\\work\\MSU\\cnv\\Results\\New_simulation_results")
Comparison<-read.table(infile,stringsAsFactors=F,sep="\t",header=FALSE)

Comparison_N<-ggplot(Comparison, aes(x=Comparison[,2], y=Comparison[,3])) +geom_point(size = 0.5)+ geom_smooth(method = "loess", se = FALSE,color="coral") + theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Comparison_F<-ggplot(Comparison, aes(x=Comparison[,2], y=Comparison[,4])) + geom_point(size = 0.5)+ geom_smooth(method = "loess", se = FALSE,color="coral") + theme_bw() +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Comparison_N<-Comparison_N+xlab("q_value") + ylab("F-measure_Nucleotide") 
Comparison_N<-Comparison_N+theme(axis.title.x = element_text(size = 10, vjust = 0), axis.title.y = element_text(size = 10, vjust = 1.5)) 
Comparison_N<-Comparison_N+theme(axis.text.x= element_text(size=5),axis.text.y= element_text(size=5))
Comparison_F<-Comparison_F+xlab("q_value") + ylab("F-measure_Region") 
Comparison_F<-Comparison_F+theme(axis.title.x = element_text(size = 10, vjust = 0), axis.title.y = element_text(size = 10, vjust = 1.5))
Comparison_F<-Comparison_F+theme(axis.text.x= element_text(size=5),axis.text.y= element_text(size=5))

pdf(paste(infile,'.pdf',sep=''),width=18,height=3)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(Comparison_N, vp = vplayout(1,1)) 
print(Comparison_F, vp = vplayout(1,2)) 
dev.off()
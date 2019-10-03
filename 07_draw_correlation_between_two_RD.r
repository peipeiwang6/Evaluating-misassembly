args = commandArgs(TRUE)
input1 = args[1] ### simulated RD
input2 = args[2] ### new RD with resampled reads
dat <- read.table(input1,head=F,sep='\t')
dat2 <- read.table(input2,head=F,sep='\t')
res <- merge(dat,dat2,by.x='V1',by.y='V1')
pdf("Correlation_between_simulated_RD_and_new_RD.pdf",width=10,height=10)
smoothScatter(res[,2],res[,3],nrpoints=0.1,xlim=c(0,5),ylim=c(0,5))
dev.off()









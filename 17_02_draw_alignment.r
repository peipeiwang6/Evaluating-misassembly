dat <- read.table('All_chr_corresponding_between_two_assemblies.txt',head=F,sep='\t',stringsAsFactor=F)
dat <- dat[dat[,9]>=97,]
dat <- dat[abs(dat[,3]-dat[,2])>=50,]
dat[,2] <- dat[,2]/1000000
dat[,3] <- dat[,3]/1000000
dat[,5] <- dat[,5]/1000000
dat[,6] <- dat[,6]/1000000
len_old <- read.table('SL2.5_merged_chr_length.txt',head=F,sep='\t',stringsAsFactor=F)
len_old[,2] <- len_old[,2]/1000000
len_new <- read.table('Sly_v4_merged_chr_length.txt',head=F,sep='\t',stringsAsFactor=F)
len_new[,2] <- len_new[,2]/1000000
pdf('All_chr_97.pdf')
a = matrix(c(0,850,0,850),nrow=2,ncol=2)
plot(a)
for(i in 1:nrow(dat)){
	if(dat[i,7]=='forward'){
		segments(dat[i,2],dat[i,5],dat[i,3],dat[i,6],col='red',lwd=0.05)
		}
	else {
		segments(dat[i,2],dat[i,5],dat[i,3],dat[i,6],col='blue',lwd=0.05)
		}
	}
old_pos = 0
for(i in 1:12) {
	old_pos = old_pos + len_old[i,2]
	abline(v=old_pos,col='grey',lty=2,lwd=0.1)
	}
new_pos = 0
for(i in 1:12) {
	new_pos = new_pos+len_new[i,2]
	abline(h=new_pos,col='grey',lty=2,lwd=0.1)
	}
dev.off()


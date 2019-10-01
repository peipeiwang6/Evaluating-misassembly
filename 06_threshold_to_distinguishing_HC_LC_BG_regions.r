'''
input 1: read depth from CNVnator
input 2: HC/LC designations from CNVnator, not that HC is refered to as duplication, LC is deletion
input 3: HC/LC designations after filtering using q0 and q-values
input 4: output1, HC with high confidence
input 5: output2, BG with high confidence
input 6: output3, LC with high confidence
'''
args = commandArgs(TRUE)
Data = args[1]
CNV = args[2]
CNV_filter = args[3]
out1 = args[4]
out2 = args[5]
out3 = args[6]

dat <- read.table(Data,head=F,stringsAsFactor=F,sep=' ')
cnv <- read.table(CNV,head=F,stringsAsFactor=F,sep='\t')
cnv_filter <- read.table(CNV_filter,head=F,stringsAsFactor=F,sep='\t')
bac <- t(t(setdiff(dat[,2],cnv[,2])))
bac <- as.data.frame(bac)
baccov <- merge(bac,dat,by.x='V1',by.y='V2')
baccov[,4] <- baccov[,4]/2  
dupcov <- cnv_filter[cnv_filter[,1]=='duplication',1:4]
delcov <- cnv_filter[cnv_filter[,1]=='deletion',1:4]
### left threshold, distringuishing LC and BG regions, to minimize the false negative rate
	F_res <- c()
	left_cut <- sort(unique(c(delcov[,4],baccov[,4])))
	left_cut <- left_cut[left_cut[]>median(delcov[,4]) & left_cut[]<median(baccov[,4])]
	for(left in left_cut){
		FN_A <- nrow(delcov[delcov[,4]>left,])/nrow(delcov)
		FN_B <- nrow(baccov[baccov[,4]<left,])/nrow(baccov)
		F_res <- rbind(F_res,c(left,FN_A,FN_B))}
	F_res <- cbind(F_res,abs(F_res[,2]-F_res[,3]))
	Left_threshold <- F_res[F_res[,4]==min(F_res[,4]),][1]
	pdf('left_threshold.pdf')
	plot(F_res[,1],F_res[,2],type='l',ylim=c(0,0.5))
	lines(F_res[,1],F_res[,3])
	dev.off()
### right threshold, distinguishing BG and HC regions
	F_res <- c()
	right_cut <- sort(unique(c(dupcov[,4],baccov[,4])))
	right_cut <- right_cut[right_cut[]>median(baccov[,4]) & right_cut[]<median(dupcov[,4])]
	for(right in right_cut){
	FN_A <- nrow(baccov[baccov[,4]>right,])/nrow(baccov)
	FN_B <- nrow(dupcov[dupcov[,4]<right,])/nrow(dupcov)
	F_res <- rbind(F_res,c(right,FN_A,FN_B))}
	F_res <- cbind(F_res,abs(F_res[,2]-F_res[,3]))
	Right_threshold <- F_res[F_res[,4]==min(F_res[,4]),][1]
	pdf('right_threshold.pdf')
	plot(F_res[,1],F_res[,2],type='l',ylim=c(0,0.5))
	lines(F_res[,1],F_res[,3])
	dev.off()
write.table(dupcov[dupcov[,4]>Right_threshold,],out1,row.names=F,col.names=F,sep='\t',quote=F)
write.table(baccov[baccov[,4]>=Left_threshold & baccov[,4]<=Right_threshold,],out2,row.names=F,col.names=F,sep='\t',quote=F)
write.table(delcov[delcov[,4]<Left_threshold,],out3,row.names=F,col.names=F,sep='\t',quote=F)

overlapping_perc <- function(left1,right1,left2,right2){
	if(as.numeric(left1-right2)*as.numeric(right1-left2) < 0){
		overlap = min(abs(left1-right2)+1, abs(right1-left2)+1,abs(left1-right1) + 1,abs(left2-right2) + 1) 
		return(c(overlap/(abs(left1-right1) + 1),overlap/(abs(left2-right2) + 1),overlap))
		}
	else {
		return(c(0,0,0))
		}
}

Region <- function(HC){
	pos1 <- regexpr(':',HC[,2])
	pos2 <- regexpr('-',HC[,2])
	HC <- cbind(HC,'Chr'=as.character(substr(HC[,2],1,pos1-1)))
	HC <- cbind(HC,'Left'=as.numeric(substr(HC[,2],pos1+1,pos2-1)))
	HC <- cbind(HC,'Right'=as.numeric(substr(HC[,2],pos2+1,nchar(HC[,2]))))
	return(HC)
	}
	
CNV_mis <- function(mis,HC){
	HC_mis <- c()
	for(i in 1:nrow(HC)){
		tem <- mis[mis[,1]==HC[i,5] & (HC[i,6]-mis[,3])*(HC[i,7]-mis[,2])<0,]
		if(nrow(tem)>=1){
			for(j in 1:nrow(tem)){
				out <- c(HC[i,1:4],as.character(HC[i,5]),HC[i,6:7],tem[j,],overlapping_perc(HC[i,6],HC[i,7],tem[j,2],tem[j,3]))
				HC_mis <- rbind(HC_mis,out)
				}
			}
		}
	HC_mis <- unique(HC_mis)
	return(HC_mis)
	}

All_2 <- read.table('All_region.txt',head=T,sep='\t',stringsAsFactor=F)
HC <- read.table('20180225_Duplication_region_with_0.7_1.76_coverage.txt',head=F,sep='\t',stringsAsFactor=F)
HC <- Region(HC)
BG <- read.table('20180225_Background_region_with_0.7_1.76_coverage.txt_modified',head=F,sep='\t',stringsAsFactor=F)
BG <- Region(BG)
LC <- read.table('20180225_Deletion_region_with_0.7_1.76_coverage.txt',head=F,sep='\t',stringsAsFactor=F)
LC <- Region(LC)
HC_mis <- CNV_mis(All_2,HC)
BG_mis <- CNV_mis(All_2,BG)
LC_mis <- CNV_mis(All_2,LC)
CNV <- rbind(HC_mis,BG_mis)
CNV <- rbind(CNV,LC_mis)
colnames(CNV) <- c('CNV','Location','Length','RD','Chr_CNV','Left','Rigth','Chr_mis_assembly','Left_mis','Right_mis','Type_mis','Dir_mis','Per_overlapping_CNV','Per_overlapping_mis','Len_overlapping')
write.table(as.matrix(as.data.frame(CNV)),'CNV_mis_dup_all.txt',row.names=F,sep='\t',quote=F)

CNV <- read.table('CNV_mis_dup_all.txt',head=T,sep='\t',stringsAsFactor=F)
CNV <- CNV[CNV[,15]>=100,]
HC_count <- aggregate(Len_overlapping ~ Type_mis, CNV[CNV[,1]=='duplication',],sum)
BG_count <- aggregate(Len_overlapping ~ Type_mis, CNV[CNV[,1]=='background',],sum)
LC_count <- aggregate(Len_overlapping ~ Type_mis, CNV[CNV[,1]=='deletion',],sum)
All_count <- aggregate(len ~ type, All_2,sum)
res <- cbind(HC_count[,1],HC_count[,2]/All_count[,2])
res <- cbind(res,BG_count[,2]/All_count[,2])
res <- cbind(res,LC_count[,2]/All_count[,2])
rownames(res) <- res[,1]
res <- res[,2:4]
out_tem <- c(1-sum(as.numeric(res[1,])),1-sum(as.numeric(res[2,])),1-sum(as.numeric(res[3,])),1-sum(as.numeric(res[4,])),1-sum(as.numeric(res[5,])),1-sum(as.numeric(res[6,])),1-sum(as.numeric(res[7,])))
colnames(res) <- c('HC','BG','LC')
res <- cbind(res,'Others'=out_tem)
res <- res[c(7,6,5,2,3,1,4),]
res2 <- cbind(HC_count[,1],HC_count[,2]/sum(HC[,3]))
res2 <- cbind(res2,BG_count[,2]/sum(BG[,3]))
res2 <- cbind(res2,LC_count[,2]/sum(LC[,3]))
row_names <- c(res2[,1],'Others')
res2 <- res2[,2:4]
out_tem <- c(1-sum(as.numeric(res2[,1])),1-sum(as.numeric(res2[,2])),1-sum(as.numeric(res2[,3])))
res2 <- rbind(res2,'Others'=out_tem)
rownames(res2) <- row_names
colnames(res2) <- c('HC','BG','LC')
res2 <- res2[c(4,1,3,2,5,6,7,8),]
pdf('Proportion_of_mis_assembly_HC_BC_LC_100bp.pdf')
par(mfrow=c(2,2))
barplot(as.numeric(res[,1]),legend=T,col=c('yellow','orange','green','blue','pink','red','purple'),las=2,beside=TRUE)
barplot(as.numeric(res[,2]),legend=T,col=c('yellow','orange','green','blue','pink','red','purple'),las=2,beside=TRUE)
barplot(as.numeric(res[,3]),legend=T,col=c('yellow','orange','green','blue','pink','red','purple'),las=2,beside=TRUE)
barplot(as.numeric(res[,4]),legend=T,col=c('yellow','orange','green','blue','pink','red','purple'),las=2,beside=TRUE)
dev.off()
pdf('Proportion_of_mis_assembly_HC_BC_LC_100bp_of_CNV.pdf')
par(mfrow=c(1,1))
barplot(res2,legend=T,col=c('yellow','orange','green','blue','pink','red','purple','grey'),ylim=c(0,1))
dev.off()


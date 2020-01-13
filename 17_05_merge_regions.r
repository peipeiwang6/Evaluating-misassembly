Link <- function(dat,List){
	dat3 <- c()
	chr <- unique(dat[,1])
	for(i in 1:length(chr)){
		dat2 <- dat[dat[,1]==chr[i],]
		pos <- c(dat2[,2],dat2[,3])
		pos <- unique(pos[order(pos)])
		for(j in 1:(length(pos)-1)){
			left <- pos[j]
			right <- pos[j+1]
			tem <- dat2[as.numeric(dat2[,3]-left) * as.numeric(dat2[,2]-right)<0,]
			if(nrow(tem)==1){
				dat3 <- rbind(dat3,c(tem[1,1],left,right,tem[1,4]))
				}
			if(nrow(tem)>1){
				Rank <- c()
				for(k in 1:nrow(tem)){
					Rank <- rbind(Rank,c(tem[k,4],which(List==as.character(tem[k,4]))))
					}
				dat3 <- rbind(dat3,c(tem[1,1],left,right,Rank[as.numeric(Rank[,2])==min(as.numeric(Rank[,2])),1]))
				}
			}
		}
	dat3 <- unique(dat3)
	return(dat3)
	}

	
All <- read.table('All_region_merged_20191205_03.txt',head=F,sep='\t',stringsAsFactor=F)
All[All[,2]<1,2] <- 1
All <- cbind(All,'Len'=All[,3]-All[,2]+1)
All_mis_dup_same_chr <- Link_subtype(All,c('mis-assembly_dup_same_chr','dup_same_chr','mis_dup_diff_chr','dup_diff_chr','mis_non_dup','Non_duplicate'))
colnames(All_mis_dup_same_chr) <- c('Chr','left','right','type','len')
write.table(as.matrix(as.data.frame(All_mis_dup_same_chr)),'All_regions.txt',row.names=F,col.names=T,sep='\t',quote=F)

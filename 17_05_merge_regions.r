Link_subtype <- function(dat,List){
	dat2 <- c()
	for(mis in List){
		tem <- dat[dat[,4]==mis,]
		tem <- unique(tem)
		tem <- tem[order(tem[,1],tem[,2]),]
		i <- 1
		while(i < nrow(tem)){
			j <- i + 1
			left <- as.numeric(tem[i,2])
			right <- as.numeric(tem[i,3])
			while(j < nrow(tem) & tem[i,1]==tem[j,1] & right-as.numeric(tem[j,2])>=0){
				right <- max(right,as.numeric(tem[j,3]))
				j <- j + 1
				}
			if(j == nrow(tem)){
				if(tem[i,1]==tem[j,1] & right-as.numeric(tem[j,2])>0){
					right <- max(right,as.numeric(tem[j,3]))
					}
				}
			dat2 <- rbind(dat2,c(tem[i,1:2],right,tem[i,4]))
			i <- j
			}
		}
	dat2 <- cbind(dat2,as.numeric(dat2[,3])-as.numeric(dat2[,2])+1)
	return(unique(dat2))
	}
	
All <- read.table('All_region_merged_20191205_03.txt',head=F,sep='\t',stringsAsFactor=F)
All[All[,2]<1,2] <- 1
All <- cbind(All,'Len'=All[,3]-All[,2]+1)
All_mis_dup_same_chr <- Link_subtype(All,c('mis-assembly_dup_same_chr','dup_same_chr','mis_dup_diff_chr','dup_diff_chr','mis_non_dup','Non_duplicate'))
colnames(All_mis_dup_same_chr) <- c('Chr','left','right','type','len')
write.table(as.matrix(as.data.frame(All_mis_dup_same_chr)),'All_regions.txt',row.names=F,col.names=T,sep='\t',quote=F)

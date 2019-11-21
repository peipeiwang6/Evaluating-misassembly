overlapping_perc <- function(left1,right1,left2,right2){
	if(as.numeric(left1-right2)*as.numeric(right1-left2) < 0){
		overlap = min(abs(left1-right2), abs(right1-left2))
		return(c(overlap/(abs(left1-right1) + 1),overlap/(abs(left2-right2) + 1),overlap))
		}
	else {
		return(c(0,0,0))
		}
}

dat1 <- read.table('All_chr_corresponding_between_two_assemblies_true_chr_pos.txt',head=F,sep='\t',stringsAsFactor=F)
subdat1 <- dat1[dat1[,1]=='NC_015438.2'|dat1[,1]=='NC_015439.2'|dat1[,1]=='NC_015440.2'|dat1[,1]=='NC_015441.2'|dat1[,1]=='NC_015442.2'|dat1[,1]=='NC_015443.2'|dat1[,1]=='NC_015444.2'|dat1[,1]=='NC_015445.2'|dat1[,1]=='NC_015446.2'|dat1[,1]=='NC_015447.2'|dat1[,1]=='NC_015448.2'|dat1[,1]=='NC_015449.2',]	
subdat1 <- subdat1[subdat1[,4]!='SL4.0ch00',]
subdat2 <- subdat1[order(subdat1[,4],subdat1[,5]),]
subdat2 <- subdat2[subdat2[,7]=='forward',]
chr_old <- unique(subdat1[,1])
chr_old <- chr_old[order(chr_old)]
chr_new <- unique(subdat1[,4])
chr_new <- chr_new[order(chr_new)]

### new corresponding regions not overlapping 
mis_dup <- c()
mis_region <- c()
for(i in 1:length(chr_old)){
	tem <- subdat1[subdat1[,1]==chr_old[i] & subdat1[,7]=='forward',]
	tem <- tem[order(tem[,2]),]
	for(j in 1:(nrow(tem)-1)){
		if(tem[j,3]-tem[j+1,2] >= 100){ ### two regions in old assembly overlap
			if(tem[j+1,4]==tem[j,4]){ ### corresponding regions in new assembly are in the same chromosome
				if(overlapping_perc(tem[j,5],tem[j,6],tem[j+1,5],tem[j+1,6])[3] < 0.5 * (tem[j,3]-tem[j+1,2])){ ### corresponding regions in new assembly can be overlapping, but no longer than overlapped region in old assembly
					b1 <- as.numeric(as.numeric(tem[j,3])*as.numeric(tem[j,5]) - as.numeric(tem[j,2])*as.numeric(tem[j,6]))/as.numeric(tem[j,3]-tem[j,2])
					b2 <- as.numeric(as.numeric(tem[j+1,3])*as.numeric(tem[j+1,5]) - as.numeric(tem[j+1,2])*as.numeric(tem[j+1,6]))/as.numeric(tem[j+1,3]-tem[j+1,2])
					a1 <- as.numeric(tem[j,6]-tem[j,5])/as.numeric(tem[j,3]-tem[j,2])
					a2 <- as.numeric(tem[j+1,6]-tem[j+1,5])/as.numeric(tem[j+1,3]-tem[j+1,2])
					region_01 <- c(tem[j+1,2],tem[j,3],tem[j+1,2]*a1+b1,tem[j,6])
					region_02 <- c(tem[j+1,2],tem[j,3],tem[j+1,5],tem[j,3]*a2+b2)
					new_01 <- subdat2[subdat2[,4]==tem[j,4] & as.numeric(subdat2[,5]-region_01[4])*as.numeric(subdat2[,6]-region_01[3])<0 & (subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]) & (subdat2[,5]!=tem[j+1,5] | subdat2[,6]!=tem[j+1,6]),]
					new_02 <- subdat2[subdat2[,4]==tem[j,4] & as.numeric(subdat2[,5]-region_02[4])*as.numeric(subdat2[,6]-region_02[3])<0 & (subdat2[,5]!=tem[j+1,5] | subdat2[,6]!=tem[j+1,6]) & (subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]),]
					if(nrow(new_01)==0 & nrow(new_02)==0){ ### corresponding regions in new assembly only have one match in old assembly
						mis_dup <- rbind(mis_dup,c(tem[j,],' '))
						mis_dup <- rbind(mis_dup,c(tem[j+1,],'mis-assembly'))
						mis_region <- rbind(mis_region,c(tem[j+1,1],region_01[1],region_02[2],'mis_assembly_tandem_duplicates'))
						}
					}
				}
			}
		}
	}
write.table(mis_dup,'Mis-assembly_by_tandem_duplication_100bp_old_half_new_soft.txt',sep='\t',quote=F,row.names=F,col.names=F)

mis_dup <- c()
for(i in 1:length(chr_old)){
	tem <- subdat1[subdat1[,1]==chr_old[i] & subdat1[,7]=='forward',]
	tem <- tem[order(tem[,2]),]
	for(j in 1:(nrow(tem)-1)){
		if(tem[j,3]-tem[j+1,2] >= 100){ ### two regions in old assembly overlap
			if(tem[j+1,4]!=tem[j,4]){ ### corresponding regions in new assembly are in the different chromosome
				if(overlapping_perc(tem[j,5],tem[j,6],tem[j+1,5],tem[j+1,6])[3] < 0.5 * (tem[j,3]-tem[j+1,2])){ ### corresponding regions in new assembly can be overlapping, but no longer than overlapped region in old assembly
					b1 <- as.numeric(as.numeric(tem[j,3])*as.numeric(tem[j,5]) - as.numeric(tem[j,2])*as.numeric(tem[j,6]))/as.numeric(tem[j,3]-tem[j,2])
					b2 <- as.numeric(as.numeric(tem[j+1,3])*as.numeric(tem[j+1,5]) - as.numeric(tem[j+1,2])*as.numeric(tem[j+1,6]))/as.numeric(tem[j+1,3]-tem[j+1,2])
					a1 <- as.numeric(tem[j,6]-tem[j,5])/as.numeric(tem[j,3]-tem[j,2])
					a2 <- as.numeric(tem[j+1,6]-tem[j+1,5])/as.numeric(tem[j+1,3]-tem[j+1,2])
					region_01 <- c(tem[j+1,2],tem[j,3],tem[j+1,2]*a1+b1,tem[j,6])
					region_02 <- c(tem[j+1,2],tem[j,3],tem[j+1,5],tem[j,3]*a2+b2)
					new_01 <- subdat2[subdat2[,4]==tem[j,4] & as.numeric(subdat2[,5]-region_01[4])*as.numeric(subdat2[,6]-region_01[3])<0 & (subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]) & (subdat2[,5]!=tem[j+1,5] | subdat2[,6]!=tem[j+1,6]),]
					new_02 <- subdat2[subdat2[,4]==tem[j,4] & as.numeric(subdat2[,5]-region_02[4])*as.numeric(subdat2[,6]-region_02[3])<0 & (subdat2[,5]!=tem[j+1,5] | subdat2[,6]!=tem[j+1,6]) & (subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]),]
					if(nrow(new_01)==0 & nrow(new_02)==0){ ### corresponding regions in new assembly only have one match in old assembly
						mis_dup <- rbind(mis_dup,c(tem[j,],' '))
						mis_dup <- rbind(mis_dup,c(tem[j+1,],'mis-assembly'))
						mis_region <- rbind(mis_region,c(tem[j+1,1],region_01[1],region_02[2],'mis_assembly_non-tandem_duplicates'))
						}
					}
				}
			}
		}
	}
write.table(mis_dup,'Mis-assembly_by_non-tandem_duplication_100bp_old.txt',sep='\t',quote=F,row.names=F,col.names=F)

mis <- c()
for(i in 1:length(chr_old)){
	tem <- subdat1[subdat1[,1]==chr_old[i] & subdat1[,7]=='forward',]
	tem <- tem[order(tem[,2],tem[,5]),]
	j <- 1
	while (j <= nrow(tem)-2){
		if(tem[j,4]==tem[j+1,4] & tem[j,4]==tem[j+2,4]){
			if(tem[j+1,5]-tem[j,6] > tem[j+2,5]-tem[j,6]){ ### there is a shift location between two assemblies
														   ### residue b = (x2*y1 - x1*y2)/(x2-x1)
				if(tem[j+2,5] > tem[j,6]){
					b0 <- as.numeric(as.numeric(tem[j,3])*as.numeric(tem[j,5]) - as.numeric(tem[j,2])*as.numeric(tem[j,6]))/as.numeric(tem[j,3]-tem[j,2])
					b1 <- as.numeric(as.numeric(tem[j+1,3])*as.numeric(tem[j+1,5]) - as.numeric(tem[j+1,2])*as.numeric(tem[j+1,6]))/as.numeric(tem[j+1,3]-tem[j+1,2])
					b2 <- as.numeric(as.numeric(tem[j+2,3])*as.numeric(tem[j+2,5]) - as.numeric(tem[j+2,2])*as.numeric(tem[j+2,6]))/as.numeric(tem[j+2,3]-tem[j+2,2])
					new_01 <- subdat2[subdat2[,4]==tem[j+1,4] & as.numeric(subdat2[,5]-tem[j+1,6])*as.numeric(subdat2[,6]-tem[j+1,5])<0 ,]
					old_01 <- subdat1[subdat1[,1]==tem[j+1,1] & as.numeric(subdat1[,2]-tem[j+1,3])*as.numeric(subdat1[,3]-tem[j+1,2])<0 ,]
					new_02 <- subdat2[subdat2[,4]==tem[j+2,4] & as.numeric(subdat2[,5]-tem[j+2,6])*as.numeric(subdat2[,6]-tem[j+2,5])<0 ,]
					old_02 <- subdat1[subdat1[,1]==tem[j+2,1] & as.numeric(subdat1[,2]-tem[j+2,3])*as.numeric(subdat1[,3]-tem[j+2,2])<0,]
					if(nrow(new_01)==1 & nrow(old_01)==1 & abs(b1-b0) > abs(b2-b0)) {
							mis <- rbind(mis,c(tem[j,],' '))
							mis <- rbind(mis,c(tem[j+1,],'mis-assembly'))
							mis <- rbind(mis,c(tem[j+2,],' '))
							mis_region <- rbind(mis_region,c(tem[j+1,1],tem[j+1,2],tem[j+1,3],'mis_assembly_non-duplicates'))
							j = j + 1
						}
					if(nrow(new_02)==1 & nrow(old_02)==1 & abs(b1-b0) < abs(b2-b0)) {
						mis <- rbind(mis,c(tem[j,],' '))
						mis <- rbind(mis,c(tem[j+1,],' '))
						mis <- rbind(mis,c(tem[j+2,],'mis-assembly'))
						mis_region <- rbind(mis_region,c(tem[j+2,1],tem[j+2,2],tem[j+2,3],'mis_assembly_non-duplicates'))
						j = j + 2
						}
					if(nrow(old_01)>1) {
						j = j + 1
						}
					if(nrow(old_02)>1) {
						j = j + 1
						}
					}
				if(tem[j+2,5] < tem[j,6]){
					j = j + 1
					}										   
				}
			}
			j = j + 1
		}
	}
write.table(mis,'Mis-assembly_100bp_old_half_new_soft.txt',sep='\t',quote=F,row.names=F,col.names=F)
write.table(mis_region,'Mis-assembly_regions.txt',sep='\t',quote=F,row.names=F,col.names=F)

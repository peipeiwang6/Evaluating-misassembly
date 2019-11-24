overlapping_perc <- function(left1,right1,left2,right2){
	if(as.numeric(left1-right2)*as.numeric(right1-left2) < 0){
		overlap = min(abs(left1-right2), abs(right1-left2))
		return(c(overlap/(abs(left1-right1) + 1),overlap/(abs(left2-right2) + 1),overlap))
		}
	else {
		return(c(0,0,0))
		}
}

setwd('D:\\Tandem_CNV\\Revision_NAR_GB')
dat <- read.table('All_chr_corresponding_between_two_assemblies_true_chr_pos.txt',head=F,sep='\t',stringsAsFactor=F)
dat_reverse <- dat[dat[,7]=='reverse',]
dat_forward <- dat[dat[,7]=='forward',]
dat1_reverse <- dat_reverse
dat1_reverse[,2] <- dat_reverse[,3]
dat1_reverse[,3] <- dat_reverse[,2]
dat1_reverse[,6] <- dat_reverse[,5]
dat1_reverse[,5] <- dat_reverse[,6]
dat1 <- rbind(dat_forward,dat1_reverse)
subdat1 <- dat1[order(dat1[,1],dat1[,2],dat1[,5]),]
subdat2 <- subdat1[order(subdat1[,4],subdat1[,5]),]
chr_old <- unique(subdat1[,1])
chr_old <- chr_old[order(chr_old)]

### new corresponding regions not overlapping 
mis_dup <- c()
mis_region <- c()
for(i in 1:length(chr_old)){
	tem <- subdat1[subdat1[,1]==chr_old[i],]
	tem <- tem[order(tem[,2]),]
	if(nrow(tem)>1) {
		for(j in 1:(nrow(tem)-1)){
			m = j+1
			while(min(tem[j,3]-tem[m,2],tem[j,3]-tem[j,2],tem[m,3]-tem[m,2]) >= 100 & m <= nrow(tem)){
				if(overlapping_perc(tem[j,5],tem[j,6],tem[m,5],tem[m,6])[3] < 0.5 * (tem[j,3]-tem[m,2])){ ### corresponding regions in new assembly can be overlapping, but no longer than overlapped region in old assembly
					b1 <- as.numeric(as.numeric(tem[j,3])*as.numeric(tem[j,5]) - as.numeric(tem[j,2])*as.numeric(tem[j,6]))/as.numeric(tem[j,3]-tem[j,2])
					b2 <- as.numeric(as.numeric(tem[m,3])*as.numeric(tem[m,5]) - as.numeric(tem[m,2])*as.numeric(tem[m,6]))/as.numeric(tem[m,3]-tem[m,2])
					a1 <- as.numeric(tem[j,6]-tem[j,5])/as.numeric(tem[j,3]-tem[j,2])
					a2 <- as.numeric(tem[m,6]-tem[m,5])/as.numeric(tem[m,3]-tem[m,2])
					x1 = max(tem[j,2],tem[m,2])
					x2 = min(tem[j,3],tem[m,3])
					region_01 <- c(x1,x2,x1*a1+b1,x2*a1+b1)
					region_02 <- c(x1,x2,x1*a2+b2,x2*a2+b2)
					new_01 <- subdat2[subdat2[,4]==tem[j,4] & as.numeric(subdat2[,5]-region_01[4])*as.numeric(subdat2[,6]-region_01[3])<0 & (subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]) & (subdat2[,5]!=tem[m,5] | subdat2[,6]!=tem[m,6]),]
					new_02 <- subdat2[subdat2[,4]==tem[j,4] & as.numeric(subdat2[,5]-region_02[4])*as.numeric(subdat2[,6]-region_02[3])<0 & (subdat2[,5]!=tem[m,5] | subdat2[,6]!=tem[m,6]) & (subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]),]
					if(nrow(new_01)==0 & nrow(new_02)==0 & min(max(tem[j,5],tem[j,6]),max(tem[m,5],tem[m,6]))-max(min(tem[j,5],tem[j,6]),min(tem[m,5],tem[m,6])) < 0.5 *min(abs(tem[j,6]-tem[j,5]),abs(tem[m,6]-tem[m,5]))){ ### corresponding regions in new assembly only have one match in old assembly
						### one region is not included in the other region in new assembly
						mis_dup <- rbind(mis_dup,c(tem[j,],' ',' '))
						if(tem[j,7]==tem[m,7]){
							if(tem[m,4]==tem[j,4]){ ### corresponding regions in new assembly are in the different chromosome
								mis_dup <- rbind(mis_dup,c(tem[m,],'mis-assembly_tandem_duplicate','forward'))
								mis_region <- rbind(mis_region,c(tem[m,1],region_01[1],region_02[2],'mis_assembly_tandem_duplicates','forward'))
								}
							else{	
								mis_dup <- rbind(mis_dup,c(tem[m,],'mis-assembly_non-tandem_duplicate','forward'))
								mis_region <- rbind(mis_region,c(tem[m,1],region_01[1],region_02[2],'mis_assembly_non-tandem_duplicates','forward'))
								}
							}
						if(tem[j,7]!=tem[m,7]){
							if(tem[m,4]==tem[j,4]){ ### corresponding regions in new assembly are in the different chromosome
								mis_dup <- rbind(mis_dup,c(tem[m,],'mis-assembly_tandem_duplicate','reverse'))
								mis_region <- rbind(mis_region,c(tem[m,1],region_01[1],region_02[2],'mis_assembly_tandem_duplicates','reverse'))
								}
							else{
								mis_dup <- rbind(mis_dup,c(tem[m,],'mis-assembly_non-tandem_duplicate','reverse'))
								mis_region <- rbind(mis_region,c(tem[m,1],region_01[1],region_02[2],'mis_assembly_non-tandem_duplicates','reverse'))
								}
							}
						}
					}
				m = m+1
				}
			}
		}
	}
write.table(mis_dup,'Mis-assembly_by_duplication_100bp_old_half_new_soft_all_direction_scaffold.txt',sep='\t',quote=F,row.names=F,col.names=F)

mis <- c()
for(i in 1:length(chr_old)){
	tem <- subdat1[subdat1[,1]==chr_old[i],]
	tem <- tem[order(tem[,2],tem[,5]),]
	if(nrow(tem)>3) {
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
								mis <- rbind(mis,c(tem[j,],' ',' '))
								mis <- rbind(mis,c(tem[j+1,],'mis-assembly','forward'))
								mis <- rbind(mis,c(tem[j+2,],' ',' '))
								mis_region <- rbind(mis_region,c(tem[j+1,1],tem[j+1,2],tem[j+1,3],'mis_assembly_non-duplicates','forward'))
								j = j + 1
							}
						if(nrow(new_02)==1 & nrow(old_02)==1 & abs(b1-b0) < abs(b2-b0)) {
							mis <- rbind(mis,c(tem[j,],' ',' '))
							mis <- rbind(mis,c(tem[j+1,],' ',' '))
							mis <- rbind(mis,c(tem[j+2,],'mis-assembly','forward'))
							mis_region <- rbind(mis_region,c(tem[j+2,1],tem[j+2,2],tem[j+2,3],'mis_assembly_non-duplicates','forward'))
							j = j + 2
							}
						if(nrow(old_01)>1) {
							j = j + 1
							}
						if(nrow(old_02)>1) {
							j = j + 1
							}
						}
					else {
						j = j + 1
						}										   
					}
				}
				j = j + 1
			}
		}
	}
write.table(mis,'Mis-assembly_100bp_old_half_new_soft_all_direction_scaffold.txt',sep='\t',quote=F,row.names=F,col.names=F)
mis_region <- unique(mis_region)
write.table(mis_region,'Mis-assembly_regions_all_direction_scaffold.txt',sep='\t',quote=F,row.names=F,col.names=F)

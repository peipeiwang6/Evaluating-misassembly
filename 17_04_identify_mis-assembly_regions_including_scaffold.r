overlapping_perc <- function(left1,right1,left2,right2){
	if(as.numeric(left1-right2)*as.numeric(right1-left2) < 0){
		overlap = min(abs(left1-right2)+1, abs(right1-left2)+1,abs(left1-right1) + 1,abs(left2-right2) + 1) 
		return(c(overlap/(abs(left1-right1) + 1),overlap/(abs(left2-right2) + 1),overlap))
		}
	else {
		return(c(0,0,0))
		}
}

Same_chr <- function(chr1,chr2,Chr){
	subchr <- Chr[Chr[,1]==chr1 & Chr[,2]==chr2,]
	if(nrow(subchr)==1) return(1)
	else return(0)
	}
	
Residue <- function(x1,x2,y1,y2){
	b = (x2*y1 - x1*y2)/(x2-x1)
	return(b)
	}	

Residue <- function(x1,x2,y1,y2){
	b = (as.numeric(x2)*as.numeric(y1) - as.numeric(x1)*as.numeric(y2))/(as.numeric(x2)-as.numeric(x1))
	return(b)
	}	

setwd('D:\\Tandem_CNV\\Revision_NAR_GB')
Chr <- read.table('Chr_in_two_assemblies.txt',head=F,sep='\t',stringsAsFactor=F)
dat <- read.table('All_chr_corresponding_between_two_assemblies_true_chr_pos_ordered.txt',head=F,sep='\t',stringsAsFactor=F)
dat <- dat[order(dat[,1],dat[,2],dat[,3],dat[,4],dat[,5]),]
### correct assembly, not duplicates
C_old <- unique(dat[,1])
NonD <- c()
Dup <- c()
for(i in 1:length(C_old)){
	tem <- dat[dat[,1]==C_old[i],]
	pos <- c(tem[,2],tem[,3])
	pos <- unique(pos[order(pos)])
	for(j in 1:(length(pos)-1)){
		left <- pos[j]
		right <- pos[j+1]
		subtem <- tem[as.numeric(tem[,3]-left)*as.numeric(tem[,2]-right) < 0,]
		if(nrow(subtem)>1){
			Dup <- rbind(Dup,c(C_old[i],left,right,'Duplicate_region'))
			}
		else {
			NonD <- rbind(NonD,c(C_old[i],left,right,'Non_duplicate'))
			}
		}
	}
write.table(Dup,'Duplication_regions.txt',row.names=F,col.names=F,quote=F,sep='\t')
write.table(NonD,'Non_duplication_regions.txt',row.names=F,col.names=F,quote=F,sep='\t')

Dup <- read.table('Duplication_regions.txt',head=F,sep='\t',stringsAsFactor=F)
NonD <- read.table('Non_duplication_regions.txt',head=F,sep='\t',stringsAsFactor=F)

subdat1 <- dat[dat[,1]=='NC_015438.2'|dat[,1]=='NC_015439.2'|dat[,1]=='NC_015440.2'|dat[,1]=='NC_015441.2'|dat[,1]=='NC_015442.2'|dat[,1]=='NC_015443.2'|dat[,1]=='NC_015444.2'|dat[,1]=='NC_015445.2'|dat[,1]=='NC_015446.2'|dat[,1]=='NC_015447.2'|dat[,1]=='NC_015448.2'|dat[,1]=='NC_015449.2',]	
subdat1 <- subdat1[subdat1[,4]!='SL4.0ch00',]
subdat1 <- subdat1[order(subdat1[,1],subdat1[,2],subdat1[,3],subdat1[,4],subdat1[,4]),]
subdat2 <- subdat1[order(subdat1[,4],subdat1[,5],subdat1[,6],subdat1[,1],subdat1[,2]),]
mis_dup_same_chr <- c()
non_mis_but_dup_same_chr <- c()
for(i in 1:nrow(Chr)){
	chr_new <- Chr[i,2]
	subdup <- Dup[Dup[,1]==Chr[i,1],]
	tem <- subdat1[subdat1[,1]==Chr[i,1] & subdat1[,4]==chr_new,]
	for(p in 1:nrow(subdup)){
		tem2 <- tem[as.numeric(tem[,3]-subdup[p,2])*as.numeric(tem[,2]-subdup[p,3])<0,]
		if(nrow(tem2)>1){
			j <- which(tem[,3]==tem2[1,3])[1]
			m <- j + 1
			while(min(tem[j,3]-tem[m,2],tem[j,3]-tem[j,2],tem[m,3]-tem[m,2]) >= 100 & m <= nrow(tem)){
				if(overlapping_perc(tem[j,5],tem[j,6],tem[m,5],tem[m,6])[3] <= 0.5 * (min(tem[j,3],tem[m,3])-max(tem[m,2],tem[j,2]))){ ### corresponding regions in new assembly can be overlapping, but no longer than overlapped region in old assembly
					b1 <- as.numeric(as.numeric(tem[j,3])*as.numeric(tem[j,5]) - as.numeric(tem[j,2])*as.numeric(tem[j,6]))/as.numeric(tem[j,3]-tem[j,2])
					b2 <- as.numeric(as.numeric(tem[m,3])*as.numeric(tem[m,5]) - as.numeric(tem[m,2])*as.numeric(tem[m,6]))/as.numeric(tem[m,3]-tem[m,2])
					a1 <- as.numeric(tem[j,6]-tem[j,5])/as.numeric(tem[j,3]-tem[j,2])
					a2 <- as.numeric(tem[m,6]-tem[m,5])/as.numeric(tem[m,3]-tem[m,2])
					x1 = max(tem[j,2],tem[m,2])
					x2 = min(tem[j,3],tem[m,3])
					region_01 <- c(x1,x2,min(x1*a1+b1,x2*a1+b1),max(x1*a1+b1,x2*a1+b1))
					region_02 <- c(x1,x2,min(x1*a2+b2,x2*a2+b2),max(x1*a2+b2,x2*a2+b2))
					new_01 <- subdat2[subdat2[,4]==tem[j,4] &  subdat2[,1] == tem[j,1] &
											as.numeric(subdat2[,5]-region_01[4])*as.numeric(subdat2[,6]-region_01[3])<0 & 
											(subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]) & 
											(subdat2[,5]!=tem[m,5] | subdat2[,6]!=tem[m,6]),]
					new_02 <- subdat2[subdat2[,4]==tem[m,4] & subdat2[,1] == tem[m,1] &
											as.numeric(subdat2[,5]-region_02[4])*as.numeric(subdat2[,6]-region_02[3])<0 & 
											(subdat2[,5]!=tem[m,5] | subdat2[,6]!=tem[m,6]) & 
											(subdat2[,5]!=tem[j,5] | subdat2[,6]!=tem[j,6]),]
					if(nrow(new_01)==0 & nrow(new_02)==0 & min(tem[j,6],tem[m,6])-max(tem[j,5],tem[m,5]) < 0.5 *abs(region_02[4]-region_02[3])){ ### corresponding regions in new assembly only have one match in old assembly
							### one region is not included in the other region in new assembly
						mis_dup_same_chr <- rbind(mis_dup_same_chr,c(tem[j,],x1,x2,'mis-assembly_dup_same_chr'))
						}
					if(nrow(new_02)>0){
						Find = 0
						new_03 <- new_02[new_02[,1]==tem[j,1],]
						if(nrow(new_03)>0){
							for(k in 1:nrow(new_03)){
								if(overlapping_perc(region_02[3],region_02[4],new_03[k,5],new_03[k,6])[3] > 0.5*abs(region_02[4]-region_02[3])){
									xk1 <- min(new_03[k,2],new_03[k,3])
									xk2 <- max(new_03[k,2],new_03[k,3])
									yk1 <- min(new_03[k,5],new_03[k,6])
									yk2 <- max(new_03[k,5],new_03[k,6])
									bk <- as.numeric(as.numeric(xk2)*as.numeric(yk1) - as.numeric(xk1)*as.numeric(yk2))/as.numeric(xk2-xk1)
									b1_2 <- as.numeric(as.numeric(max(tem[j,2],tem[j,3]))*as.numeric(min(tem[j,5],tem[j,6])) - as.numeric(min(tem[j,2],tem[j,3]))*as.numeric(max(tem[j,5],tem[j,6])))/as.numeric(abs(tem[j,2]-tem[j,3]))
									b2_2 <- as.numeric(as.numeric(max(tem[m,2],tem[m,3]))*as.numeric(min(tem[m,5],tem[m,6])) - as.numeric(min(tem[m,2],tem[m,3]))*as.numeric(max(tem[m,5],tem[m,6])))/as.numeric(abs(tem[m,2]-tem[m,3]))
									if(abs(bk - b1_2) < abs(b2_2 - b1_2)){
										Find = 1
										}
									}
								}
							}
						if(Find == 0) {
							mis_dup_same_chr <- rbind(mis_dup_same_chr,c(tem[j,],x1,x2,'mis-assembly_dup_same_chr'))
							}
						if(Find == 1) {
							non_mis_but_dup_same_chr <- rbind(non_mis_but_dup_same_chr,c(tem[j,],x1,x2,'dup_same_chr'))
							}
						}
					}
				if(overlapping_perc(tem[j,5],tem[j,6],tem[m,5],tem[m,6])[3] > 0.5 * (min(tem[j,3],tem[m,3])-max(tem[m,2],tem[j,2]))){
					x1 = max(tem[j,2],tem[m,2])
					x2 = min(tem[j,3],tem[m,3])
					non_mis_but_dup_same_chr <- rbind(non_mis_but_dup_same_chr,c(tem[j,],x1,x2,'dup_same_chr'))
					}
				m = m+1
				}
			}
		}
	}
mis_dup_same_chr <- unique(mis_dup_same_chr)
non_mis_but_dup_same_chr <- unique(non_mis_but_dup_same_chr)
write.table(mis_dup_same_chr,'Mis-assembly_dup_same_chr.txt',sep='\t',quote=F,row.names=F,col.names=F)
write.table(non_mis_but_dup_same_chr,'Non-misassembly_Dup_same_chr.txt',sep='\t',quote=F,row.names=F,col.names=F)

mis_non_dup <- c()  ### not solved yet, have some issue
for(i in 1:length(C_old)){
	subdup <- NonD[NonD[,1]==C_old[i],]
	subdat <- dat[dat[,1]==C_old[i],]
	region <- rbind(Dup[Dup[,1]==C_old[i],],NonD[NonD[,1]==C_old[i],])
	region <- region[order(region[,2]),]
	for(p in 1:nrow(region)){
		tem <- subdat[as.numeric(subdat[,3]-region[p,2])*as.numeric(subdat[,2]-region[p,3])<0,]
		if(nrow(tem)==1){
			checked = 0
			if(nrow(Chr[Chr[,1]==tem[1,1],])==1 & nrow(Chr[Chr[,2]==tem[1,4],])==1){  
				if(Same_chr(tem[1,1],tem[1,4],Chr)==0){ ### old and new are not the same chr
					if(overlapping_perc(region[p,2],region[p,3],tem[1,2],tem[1,3])[3] > 0.5*(tem[1,3]-tem[1,2]) & region[p,4]=='Non_duplicate'){
						mis_non_dup <- rbind(mis_non_dup,c(tem[1,],'mis_non_dup'))
						}
					else {
						checked = 1
						}
					}
				}
			if(checked == 0){
				for(o in 1:nrow(tem)){
					if(Same_chr(tem[o,1],tem[o,4],Chr)==1){
						subtem <- subdat[subdat[,4]==tem[o,4],]
						subtem <- subtem[order(subtem[,2]),]
						m <- which(subtem[,2]==tem[o,2])[1]
						e <- which(subtem[,3]==max(subtem[1:(m-1),3]))
						if(length(e)>1){
							residue <- c()
							for(x in 1:length(e)){
								residue <- rbind(residue,c(e[x],abs(Residue(subtem[e[x],2],subtem[e[x],3],subtem[e[x],5],subtem[e[x],6]))))
								}
							e <- residue[residue[,2]==min(residue[,2]),1]
							}
						f <- which(subtem[,2]==min(subtem[(m+1):nrow(subtem),2]))
						if(length(f)>1){
							residue <- c()
							for(x in 1:length(f)){
								residue <- rbind(residue,c(f[x],abs(Residue(subtem[f[x],2],subtem[f[x],3],subtem[f[x],5],subtem[f[x],6]))))
								}
							f <- residue[residue[,2]==min(residue[,2]),1]
							}
						if(length(f) > 0){
							if(m >= 2 &  m <= nrow(subtem)-1 &  e >= 1 & (subtem[m,7]=='forward' | subtem[f,7]=='forward')) {
								if(subtem[m,4]==subtem[e,4] & subtem[e,4]==subtem[f,4] & subtem[m,5] > subtem[f,6] & subtem[m,3] < subtem[f,2]) { ### there is a shift location between two assemblies
									if(subtem[f,5] > subtem[e,6]){
										b0 <- as.numeric(as.numeric(subtem[e,3])*as.numeric(subtem[e,5]) - as.numeric(subtem[e,2])*as.numeric(subtem[e,6]))/as.numeric(subtem[e,3]-subtem[e,2])
										b1 <- as.numeric(as.numeric(subtem[m,3])*as.numeric(subtem[m,5]) - as.numeric(subtem[m,2])*as.numeric(subtem[m,6]))/as.numeric(subtem[m,3]-subtem[m,2])
										b2 <- as.numeric(as.numeric(subtem[f,3])*as.numeric(subtem[f,5]) - as.numeric(subtem[f,2])*as.numeric(subtem[f,6]))/as.numeric(subtem[f,3]-subtem[f,2])
										new_01 <- subtem[subtem[,4]==subtem[e,4] & as.numeric(subtem[,5]-subtem[m,6])*as.numeric(subtem[,6]-subtem[m,5])<0 ,]
										new_02 <- subtem[subtem[,4]==subtem[e,4] & as.numeric(subtem[,5]-subtem[f,6])*as.numeric(subtem[,6]-subtem[f,5])<0 ,]
										old_01 <- subtem[as.numeric(subtem[,2]-subtem[m,3])*as.numeric(subtem[,3]-subtem[m,2])<0 ,]
										old_02 <- subtem[subtem[,4]==subtem[f,4] & as.numeric(subtem[,2]-subtem[f,3])*as.numeric(subtem[,3]-subtem[f,2])<0,]
										if(nrow(new_01[abs(new_01[,8])>=0.5*abs(subtem[m,8]),])==1 & nrow(old_01[abs(old_01[,8])>=0.5*abs(subtem[m,8]),])==1 & abs(b1-b0) > abs(b2-b0)) {
											mis_non_dup <- rbind(mis_non_dup,c(subtem[m,],'mis_non_dup'))
											checked = 1
											}
										}
									}
								}
							}
						if(checked == 0){
							subtem <- subdat[subdat[,4]==tem[o,4],]
							subtem <- subtem[order(subtem[,2]),]
							m <- which(subtem[,2]==tem[o,2])[1]
							f <- which(subtem[,3]==max(subtem[1:(m-1),3]))[1]
							if(length(f)>1){
								residue <- c()
								for(x in 1:length(f)){
									residue <- rbind(residue,c(f[x],abs(Residue(subtem[f[x],2],subtem[f[x],3],subtem[f[x],5],subtem[f[x],6]))))
									}
								f <- residue[residue[,2]==min(residue[,2]),1]
								}
							e <- f-1
							if(length(f) > 0){
								if(m >= 3 & m <= nrow(subtem) & e >= 1 & (subtem[m,7]=='forward' | subtem[f,7]=='forward')){
									if(subtem[e,4]==subtem[f,4] & subtem[e,4]==subtem[m,4] & subtem[f,5] > subtem[m,6] & subtem[f,3] < subtem[m,2]) { ### there is a shift location between two assemblies
										if(subtem[m,5] > subtem[e,6]){
											b0 <- as.numeric(as.numeric(subtem[e,3])*as.numeric(subtem[e,5]) - as.numeric(subtem[e,2])*as.numeric(subtem[e,6]))/as.numeric(subtem[e,3]-subtem[e,2])
											b1 <- as.numeric(as.numeric(subtem[f,3])*as.numeric(subtem[f,5]) - as.numeric(subtem[f,2])*as.numeric(subtem[f,6]))/as.numeric(subtem[f,3]-subtem[f,2])
											b2 <- as.numeric(as.numeric(subtem[m,3])*as.numeric(subtem[m,5]) - as.numeric(subtem[m,2])*as.numeric(subtem[m,6]))/as.numeric(subtem[m,3]-subtem[m,2])
											new_01 <- subtem[subtem[,4]==subtem[e,4] & as.numeric(subtem[,5]-subtem[f,6])*as.numeric(subtem[,6]-subtem[f,5])<0 ,]
											new_02 <- subtem[subtem[,4]==subtem[e,4] & as.numeric(subtem[,5]-subtem[m,6])*as.numeric(subtem[,6]-subtem[m,5])<0 ,]
											old_01 <- subtem[as.numeric(subtem[,2]-subtem[f,3])*as.numeric(subtem[,3]-subtem[f,2])<0 ,]
											old_02 <- subtem[subtem[,4]==subtem[m,4] & as.numeric(subtem[,2]-subtem[m,3])*as.numeric(subtem[,3]-subtem[m,2])<0,]
											if(nrow(new_02[abs(new_02[,8])>=0.5*abs(subtem[m,8]),])==1 & nrow(old_02[abs(old_02[,8])>=0.5*abs(subtem[m,8]),])==1 & abs(b1-b0) < abs(b2-b0)) {
												mis_non_dup <- rbind(mis_non_dup,c(subtem[m,],'mis_non_dup'))
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
mis_non_dup <- unique(mis_non_dup)	
write.table(mis_non_dup,'Mis-assembly_but_non_duplicates.txt',row.names=F,col.names=F,quote=F,sep='\t')	

mis_dup_diff_chr <- c()
dup_diff_chr <- c()
subdat3 <- dat[order(dat[,4],dat[,5],dat[,6],dat[,1],dat[,2]),]
for(i in 1:length(C_old)){
	subdup <- Dup[Dup[,1]==C_old[i],]
	subdup <- subdup[order(subdup[,2]),]
	subdat <- dat[dat[,1]==C_old[i],]
	for(p in 1:nrow(subdup)){
		tem <- subdat[as.numeric(subdat[,3]-subdup[p,2])*as.numeric(subdat[,2]-subdup[p,3])<0,]
		if(nrow(tem)>1){
			for(j in 1:(nrow(tem)-1)){
				m <- j + 1
				while(m <= nrow(tem) & min(tem[j,3]-tem[m,2],tem[j,3]-tem[j,2],tem[m,3]-tem[m,2]) >= 100){
					if(tem[j,4] != tem[m,4]){ ### corresponding regions in new assembly can be overlapping, but no longer than overlapped region in old assembly
						b1 <- as.numeric(as.numeric(tem[j,3])*as.numeric(tem[j,5]) - as.numeric(tem[j,2])*as.numeric(tem[j,6]))/as.numeric(tem[j,3]-tem[j,2])
						b2 <- as.numeric(as.numeric(tem[m,3])*as.numeric(tem[m,5]) - as.numeric(tem[m,2])*as.numeric(tem[m,6]))/as.numeric(tem[m,3]-tem[m,2])
						a1 <- as.numeric(tem[j,6]-tem[j,5])/as.numeric(tem[j,3]-tem[j,2])
						a2 <- as.numeric(tem[m,6]-tem[m,5])/as.numeric(tem[m,3]-tem[m,2])
						x1 = max(tem[j,2],tem[m,2])
						x2 = min(tem[j,3],tem[m,3])
						region_01 <- c(x1,x2,min(x1*a1+b1,x2*a1+b1),max(x1*a1+b1,x2*a1+b1))
						region_02 <- c(x1,x2,min(x1*a2+b2,x2*a2+b2),max(x1*a2+b2,x2*a2+b2))
						new_01 <- subdat3[subdat3[,1] == tem[j,1] &
												as.numeric(subdat3[,5]-region_01[4])*as.numeric(subdat3[,6]-region_01[3])<0 & 
												(subdat3[,5]!=tem[j,5] | subdat3[,6]!=tem[j,6]) & 
												(subdat3[,5]!=tem[m,5] | subdat3[,6]!=tem[m,6]),]
						new_02 <- subdat3[subdat3[,4]==tem[m,4] &
												as.numeric(subdat3[,5]-region_02[4])*as.numeric(subdat3[,6]-region_02[3])<0 & 
												(subdat3[,5]!=tem[m,5] | subdat3[,6]!=tem[m,6]) & 
												(subdat3[,5]!=tem[j,5] | subdat3[,6]!=tem[j,6]),]
						if(nrow(new_01)==0 & nrow(new_02)==0){ ### corresponding regions in new assembly only have one match in old assembly
								### one region is not included in the other region in new assembly
							mis_dup_diff_chr <- rbind(mis_dup_diff_chr,c(tem[j,],tem[m,4],x1,x2,'mis_dup_diff_chr'))
							}
						if(nrow(new_02)>0){
							Find = 0
							if(nrow(new_02)>0){
								for(k in 1:nrow(new_02)){
									if(overlapping_perc(region_02[3],region_02[4],new_02[k,5],new_02[k,6])[3] > 0.5*abs(region_02[4]-region_02[3])){
										Find = 1
										}
									}
								}
							if(Find == 0) {
								mis_dup_diff_chr <- rbind(mis_dup_diff_chr,c(tem[j,],tem[m,4],x1,x2,'mis_dup_diff_chr'))
								}
							if(Find == 1) {
								dup_diff_chr <- rbind(dup_diff_chr,c(tem[j,],tem[m,4],x1,x2,'dup_diff_chr'))
								}
							}
						}
					m = m+1
					}
				}
			}
		}
	}
dup_diff_chr <- unique(dup_diff_chr)
mis_dup_diff_chr <- unique(mis_dup_diff_chr)
write.table(dup_diff_chr,'Dup_in_diff_chr.txt',row.names=F,col.names=F,quote=F,sep='\t')	
write.table(mis_dup_diff_chr,'Mis_assembly_Dup_in_diff_chr.txt',row.names=F,col.names=F,quote=F,sep='\t')	

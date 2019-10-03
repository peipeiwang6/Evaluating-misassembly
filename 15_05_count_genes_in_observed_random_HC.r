args=commandArgs(TRUE)
observed_HC <- args[1] # '20180225_Duplication_region_with_0.7_1.76_coverage.txt'
path <- args[2] # path where the random selected HC regions were saved
tandem_file <- args[3] # sl.tandem.genes_list
proximal_file <- args[4] # sl.proximal.genes_list
sm_file <- args[5] # specialized metabolism gene list, Sol_SM_genes_Beth.csv
output_file <- args[6] # specialized metabolism gene list, Sol_SM_genes_Beth.csv

random = paste(path,'Randomly_HC_region_',sep='')
tandem <- read.table(tandem_file,head=F,stringsAsFactors=F)
pro <- read.table(proximal_file,head=F,stringsAsFactors=F)
tandem <- rbind(tandem,pro)
sm <- read.csv(sm_file,head=T,stringsAsFactors=F)
res <- c()

dat_gene <- read.table(paste(observed_HC,'_duplication_genes_perbp.txt',sep=''),head=F,sep='\t',stringsAsFactors=F)
dat_pseu <- read.table(paste(observed_HC,'_duplication_pseudogenes_perbp.txt',sep=''),head=F,sep='\t',stringsAsFactors=F)
dat_repeat <- read.table(paste(observed_HC,'_duplication_repeats_perbp.txt',sep=''),head=F,sep='\t',stringsAsFactors=F)
ngene <- nrow(dat_gene)
ngene50 <- nrow(dat_gene[dat_gene[,7]>=0.5,])
nprotein50 <- nrow(dat_gene[(dat_gene[,3]=='mRNA') & dat_gene[,7]>=0.5,])
protein50 <- dat_gene[(dat_gene[,3]=='mRNA') & dat_gene[,7]>=0.5,]
nncRNA50 <- nrow(dat_gene[dat_gene[,3]=='ncRNA' & dat_gene[,7]>=0.5,])
ntRNA50 <- nrow(dat_gene[dat_gene[,3]=='tRNA' & dat_gene[,7]>=0.5,])
nmisc_RNA50 <- nrow(dat_gene[dat_gene[,3]=='transcript' & dat_gene[,7]>=0.5,])
nmicroRNA50 <- nrow(dat_gene[dat_gene[,3]=='primary_transcript' & dat_gene[,7]>=0.5,])
npseu50 <- nrow(dat_pseu[dat_pseu[,7]>=0.5,])
nrepeat50 <- nrow(unique(dat_repeat[dat_repeat[,7]>=0.5,4:6]))
ntandem50 <- length(intersect(protein50[,2],tandem[,1]))
nsm50 <- length(intersect(protein50[,2],sm[,1]))
nsmtan50 <- length(intersect(tandem[,1],intersect(protein50[,2],sm[,1])))

res <- rbind(res,c('Observed',ngene50,nprotein50,nncRNA50,ntRNA50,nmisc_RNA50,nmicroRNA50,npseu50,nrepeat50,ntandem50,nsm50,nsmtan50))

for(i in 1:10000){
tryCatch( 
{dat_pseu <- read.table(paste(random,i,'.txt_duplication_pseudogenes_perbp.txt',sep=''),head=F,sep='\t',stringsAsFactors=F)
	npseu <- nrow(dat_pseu)
	npseu50 <- nrow(dat_pseu[dat_pseu[,7]>=0.5,])},
	error = function(e) {print(paste('Errors with ',random,i,'.txt_duplication_pseudogenes_perbp.txt',sep=''));NaN})

tryCatch( 
{dat_repeat <- read.table(paste(random,i,'.txt_duplication_repeats_perbp.txt',sep=''),head=F,sep='\t',stringsAsFactors=F)
	nrepeat <- nrow(dat_repeat)
	nrepeat50 <- nrow(dat_repeat[dat_repeat[,7]>=0.5,])},
	error = function(e) {print(paste('Errors with ',random,i,'.txt_duplication_repeats_perbp.txt',sep=''));NaN})

tryCatch(   
{
	tryCatch( 
		{dat_gene <- read.table(paste(random,i,'.txt_duplication_genes_perbp.txt',sep=''),head=F,sep='\t',stringsAsFactors=F)
			ngene <- nrow(dat_gene)
			ngene50 <- nrow(dat_gene[dat_gene[,7]>=0.5,])
			nprotein50 <- nrow(dat_gene[(dat_gene[,3]=='mRNA') & dat_gene[,7]>=0.5,])
			protein50 <- dat_gene[(dat_gene[,3]=='mRNA') & dat_gene[,7]>=0.5,]
			nncRNA50 <- nrow(dat_gene[dat_gene[,3]=='ncRNA' & dat_gene[,7]>=0.5,])
			ntRNA50 <- nrow(dat_gene[dat_gene[,3]=='tRNA' & dat_gene[,7]>=0.5,])
			nmisc_RNA50 <- nrow(dat_gene[dat_gene[,3]=='transcript' & dat_gene[,7]>=0.5,])
			nmicroRNA50 <- nrow(dat_gene[dat_gene[,3]=='primary_transcript' & dat_gene[,7]>=0.5,])
			npseu50 <- nrow(dat_pseu[dat_pseu[,7]>=0.5,])
			nrepeat50 <- nrow(dat_repeat[dat_repeat[,7]>=0.5,])
			ntandem50 <- length(intersect(protein50[,2],tandem[,1]))
			nsm50 <- length(intersect(protein50[,2],sm[,1]))
			nsmtan50 <- length(intersect(tandem[,1],intersect(protein50[,2],sm[,1])))
			res <- rbind(res,c('Randem',ngene50,nprotein50,nncRNA50,ntRNA50,nmisc_RNA50,nmicroRNA50,npseu50,nrepeat50,ntandem50,nsm50,nsmtan50))},
		error = function(e) {print(paste('Errors with ',random,i,'.txt_duplication_genes_perbp.txt',sep=''));NaN})
},
error = function(e) {
res <- rbind(res,c(paste('Randem',i,sep=''),0,0,0,0,0,0,0,npseu,npseu50,nrepeat,nrepeat50,0,0,0,0,0,0,0,0))})
}
colnames(res) <- c('Result','gene','protein-coding50','ncRNA50','tRNA50','misc_RNA50','microRNA50','pseu50','repeat50','tandem50','sm50','smtan50')
write.table(res,output_file,quote=F,sep='\t',row.names=F)
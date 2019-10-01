import sys,os
gene = open('Solanum_lycopersicum_GCF_000188115.3_S.lycopersicum.2.50_genomic.gff','r').readlines() ### GFF file from tomato genome resouces of NCBI v2.5
repe = open('ITAG2.4_repeats.gff3','r').readlines() ### annotation of transposon element, downloaded from Solanaceae Genomics Network ITAG2.4 release 
chr = open("ITAG2.4_chromosome_length.txt",'r').readlines() ### length of chromosomes
seg = open("Sly_duplication_deletion_100bp_0.7_1.76.txt",'r').readlines() ### HC regions with high confidence
dic = open("Directionary_chromosome.txt",'r').readlines() ### correlation of chromosome names between genome assembly version NCBI v2.5 and Solanaceae Genomics Network ITAG2.4 release 
tan = open('Sly_tandem_proximal_genes.txt','r').readlines()  ### tandem and proximal gene list identified using MCScanX-transposed (Wang et al 2013. Bioinformatics)
gff = open('sl.gff','r').readlines()  ### simplied GFF file
out = open('Sly_duplication_deletion_100bp_0.7_1.76_topography.txt','w')

D = {}
for inl in dic:
	D[inl.split('\t')[0]] = inl.split('\t')[1].strip()

L = {}  ### for chromosome length
for inl in chr:
	L[inl.split()[0]] = int(inl.split()[2])

for inl in gene:
	if not inl.startswith('#') and inl.split('\t')[2]=='gene' and inl.split('\t')[0] in D:
		ch = D[inl.split('\t')[0]]
		L[ch] = max(int(inl.split('\t')[4]),L[ch])

for inl in repe:
	if not inl.startswith('SL2.50ch00') and not inl.startswith('#'):
		ch = inl.split('\t')[0]
		L[ch] = max(int(inl.split('\t')[4]),L[ch])

Res = {}
for inl in chr:
	if inl.split()[0] not in Res:
		Res[inl.split()[0]] = {}
	for i in range(1,L[inl.split()[0]]+1):
		Res[inl.split()[0]][i] = [0,0,0,0,0,0]
		if i%10000==0:
			print inl.split()[0] + ': %s'%(i)

for inl in gene:
	if not inl.startswith('#') and inl.split('\t')[2]=='gene' and inl.split('\t')[0] in D:
		ch = D[inl.split('\t')[0]]
		left = int(inl.split('\t')[3])
		right = int(inl.split('\t')[4])
		for i in range(left,right+1):
			Res[ch][i][0] = 1

GFF = {}
for inl in gff:
	if inl.split('\t')[0] in D:
		GFF[inl.split('\t')[1]] = [D[inl.split('\t')[0]],int(inl.split('\t')[2]),int(inl.split('\t')[3].strip())]

T = {}  ### location for Tandem and proximal genes
for inl in tan:
	if inl.strip() in GFF:
		ch = GFF[inl.strip()][0]
		left = GFF[inl.strip()][1]
		right = GFF[inl.strip()][2]
		T[inl.strip()] = 1
		for i in range(left,right+1):
			Res[ch][i][1] = 1

for g in GFF.keys():
	if g not in T:
		ch = GFF[g][0]
		left = GFF[g][1]
		right = GFF[g][2]
		T[inl.strip()] = 1
		for i in range(left,right+1):
			Res[ch][i][2] = 1

for inl in repe:
	if not inl.startswith('SL2.50ch00') and not inl.startswith('#'):
		ch = inl.split('\t')[0]
		left = int(inl.split('\t')[3])
		right = int(inl.split('\t')[4])
		for i in range(left,right+1):
			Res[ch][i][3] = 1

for inl in seg:
	if inl.split('\t')[1].split(':')[0] in D:
		ch = D[inl.split('\t')[1].split(':')[0]]
		left = int(inl.split('\t')[1].split(':')[1].split('-')[0])
		right = int(inl.split('\t')[1].split(':')[1].split('-')[1])
		if inl.startswith('duplication'):
			for i in range(left,right+1):
				Res[ch][i][4] = 1
		if inl.startswith('deletion'):
			for i in range(left,right+1):
				Res[ch][i][5] = 1

out.write("chr\tnumber\tgene_density\ttandem_density\tnon-tandem_density\trepeat_density\thigh-coverage_density\tlow-coverage_density\n")
for inl in chr:
	ch = inl.split()[0]
	leng = L[ch]
	for i in range(1,leng/500000+1):
		n1 = 0
		n2 = 0
		n3 = 0
		n4 = 0
		n5 = 0
		n6 = 0
		for j in range((i-1)*500000+1,i*500000+1):
			n1 += Res[ch][j][0]
			n2 += Res[ch][j][1]
			n3 += Res[ch][j][2]
			n4 += Res[ch][j][3]
			n5 += Res[ch][j][4]
			n6 += Res[ch][j][5]
		out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(ch,i,float(n1)/float(500000),float(n2)/float(500000),float(n3)/float(500000),float(n4)/float(500000),float(n5)/float(500000),float(n6)/float(500000)))
		out.flush()

out.close()
















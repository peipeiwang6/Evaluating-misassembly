import sys,os

genome_file = sys.argv[1] #Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fna.longest.mod.fa
gff_file = sys.argv[2] #Solanum_lycopersicum_GCF_000188115.3_S.lycopersicum.2.50_genomic.gff
repeats_gff_file =sys.argv[3] #ITAG2.4_repeats.gff3
repeats_gff_file_2 = sys.argv[4] #Repeat_sequences_2.5genome_all_best_mid_final.gff, repeats in scaffold sequences in NCBI version 2.5
dic_file = sys.argv[5] #Dictionary of chromosome names between NCBI and SGN versions
tandem_file = sys.argv[6] #tandem gene file
proximal_file = sys.argv[7] #proximal gene file

seq = open(genome_file,'r').readlines()
S = {}
L = {}  ### for chromosome length
x = 0
while x < len(seq)/2:
	S[seq[2*x].strip()[5:]] = seq[2*x+1].strip()
	L[seq[2*x].strip()[5:]] = len(seq[2*x+1].strip())
	x += 1

def LWN(seqe):
	n = 0
	for s in seqe:
		if s !='N' and s != 'n':
			n += 1
	return n


gff = open(gff_file,'r').readlines()
repe = open(repeats_gff_file,'r').readlines()
repe_2 = open(repeats_gff_file_2,'r').readlines()
dic = open(dic_file,'r').readlines()
tandem = open(tandem_file,'r').readlines()
proximal = open(proximal_file,'r').readlines()
D = {}
for inl in dic:
	D[inl.split('\t')[0]] = inl.split('\t')[1].strip()

CDS = {} ### CDS[cds] = mRNA
M = {}  ### M[mRNA] = gene
G_pre = {}  ### for gene, Note that some gene would be overlapping with each others
x = 0
while x < len(gff):
	if not gff[x].startswith('#'):
		if gff[x].split('\t')[2] == 'gene' and 'pseudo=true' not in gff[x]:
			gene = gff[x].split('\t')[8].split(';')[0].split('ID=')[1]
			ch = gff[x].split('\t')[0]
			left = int(gff[x].split('\t')[3])
			right = int(gff[x].split('\t')[4])
			if ch != 'NC_007898.3': # this is chloroplast chromosome
				if ch not in G_pre:
					G_pre[ch] = {}
				G_pre[ch][left] = [gene,right]
		if 'mRNA' == gff[x].split('\t')[2]:
			tem = gff[x].split('\t')[8].strip()
			M[tem.split('ID=')[1].split(';')[0]] = tem.split('Parent=')[1].split(';')[0]
		if 'CDS' == gff[x].split('\t')[2]:
			tem = gff[x].split('\t')[8].strip()
			CDS[tem.split('protein_id=')[1].split(';')[0]] = tem.split('Parent=')[1].split(';')[0]
	x += 1

nn = 0
while nn < 50:
	nn += 1
	for ch in sorted(G_pre.keys()):
		left = sorted(G_pre[ch].keys())
		x = 0
		while x < len(left)-1:
			y = x + 1
			while y < len(left) and G_pre[ch][left[x]][1] >= left[y]:
				G_pre[ch][left[x]][1] = G_pre[ch][left[y]][1]
				G_pre[ch].pop(left[y])
				left.pop(y)
				y += 1
			x = y

G = {} ### for gene
for ch in sorted(G_pre.keys()):
	if ch not in G:
		G[ch] = {}
	for left in sorted(G_pre[ch].keys()):
		G[ch][G_pre[ch][left][0]] = [left,G_pre[ch][left][1]]

chromosome = open(dic_file,'r').readlines()
C = {}
for inl in chromosome:
	C[inl.strip().split('\t')[1]] = inl.split('\t')[0]

R_pre = {}  ### for repeat, Note that some repeat would be overlapping with each others
x = 0
while x < len(repe): ### repeat in assembled chromosomes
	if not repe[x].startswith('SL2.50ch00') and not repe[x].startswith('#'): 
		gene = repe[x].split('\t')[2] + '_%s'%(x) 
		ch = C[repe[x].split('\t')[0]]
		if ch not in R_pre:
			R_pre[ch] = {}
		left = int(repe[x].split('\t')[3])
		right = int(repe[x].split('\t')[4].strip())
		R_pre[ch][left] = [gene,right]
	x += 1

y = 0
while y < len(repe_2):   ### repeat in scaffolds
	inl = repe_2[y]
	tem = inl.split('\t')
	gene = tem[2] + '_%s'%(x+y) 
	ch = tem[0]
	if ch not in R_pre:
		R_pre[ch] = {}
	left = int(tem[3])
	right = int(tem[4].strip())
	R_pre[ch][left] = [gene,right]
	y += 1

### merge overlapping repeat
nn = 0
total_nn = 0
for ch in sorted(R_pre.keys()):
	for left in sorted(R_pre[ch].keys()):
		total_nn += len(R_pre[ch][left])

while total_nn != nn:
	nn = total_nn
	for ch in sorted(R_pre.keys()):
		left = sorted(R_pre[ch].keys())
		x = 0
		while x < len(left)-1:
			y = x + 1
			while y < len(left) and R_pre[ch][left[x]][1] >= left[y]:
				R_pre[ch][left[x]][1] = R_pre[ch][left[y]][1]
				R_pre[ch].pop(left[y])
				left.pop(y)
				y += 1
			x = y
	total_nn = 0
	for ch in sorted(R_pre.keys()):
		for left in sorted(R_pre[ch].keys()):
			total_nn += len(R_pre[ch][left])


R = {} ### for repeat
for ch in sorted(R_pre.keys()):
	if ch not in R:
		R[ch] = {}
	for left in sorted(R_pre[ch].keys()):
		R[ch][R_pre[ch][left][0]] = [left,R_pre[ch][left][1]]

T = {} ### for tandem
for inl in tandem:
	if not inl.startswith('Duplicate'):
		T[M[CDS[inl.split('\t')[0]]]] = 1

for inl in proximal:  ### merge proximal with tandem
	if not inl.startswith('Duplicate'):
		T[M[CDS[inl.split('\t')[0]]]] = 1

def GC_content(seq):
	seq = seq.upper()
	c_num = seq.count('C')
	g_num = seq.count('G')
	a_num = seq.count('A')
	t_num = seq.count('T')
	if c_num+g_num+a_num+t_num!=0:
		return float(c_num+g_num)/(c_num+g_num+a_num+t_num)
	else:
		return 'NaN'

dup = open(path + sys.argv[1],'r').readlines()
type = sys.argv[2]
out = open('GC_content_and_densities_exclude_N','w')
out.write('type\tregion\tWithin_or_Neighbor\tGC\tGene density\tRepeat density\tTandem density\n')
for inl in dup:
	tem = inl.split()
	ch = tem[1].split(':')[0]
	region_left = int(tem[1].split(':')[1].split('-')[0])-1   ### region left
	region_right = int(tem[1].split(':')[1].split('-')[1])   ### region right
	if region_right > len(S[ch]):
		region_right = len(S[ch])
	seq = S[ch][region_left:region_right]  ### region sequence
	GC = GC_content(seq)
	gene_length = 0
	repeat_length = 0
	tandem_length = 0
	region_length = LWN(S[ch][region_left:region_right])
	if region_length == 0:
		out.write('%s\t%s:%s_%s\tWithin\t%s\t%s\t%s\t%s\n'%(type,ch,region_left,region_right,GC,0,0,0))
	else:
		### get total gene length for region sequences
		if ch in G:
			for gene in G[ch].keys():
				if (region_left-G[ch][gene][1])*(region_right-G[ch][gene][0]) < 0 or region_left-G[ch][gene][0]==0 or region_right-G[ch][gene][1]==0:
					length = LWN(S[ch][G[ch][gene][0]:G[ch][gene][1]])
					distance = min(LWN(S[ch][min(region_left,G[ch][gene][1]):max(region_left,G[ch][gene][1])]),LWN(S[ch][min(region_right,G[ch][gene][0]):max(region_right,G[ch][gene][0])]),length,region_length)
					gene_length += distance
					if gene in T:
						tandem_length += distance
		### get total repeat length for region sequences
		if ch in R:
			for gene in R[ch].keys():
				if (region_left-R[ch][gene][1])*(region_right-R[ch][gene][0]) < 0 or region_left-R[ch][gene][0]==0 or region_right-R[ch][gene][1]==0:
					length = LWN(S[ch][R[ch][gene][0]:R[ch][gene][1]])
					distance = min(LWN(S[ch][min(region_left,R[ch][gene][1]):max(region_left,R[ch][gene][1])]),LWN(S[ch][min(region_right,R[ch][gene][0]):max(region_right,R[ch][gene][0])]),length,region_length)
					repeat_length += distance
		out.write('%s\t%s:%s_%s\tWithin\t%s\t%s\t%s\t%s\n'%(type,ch,region_left,region_right,GC,float(gene_length)/region_length,float(repeat_length)/region_length,float(tandem_length)/region_length))
		out.flush()

for l in [500,1000,2000,4000,8000,16000,32000]:
	for inl in dup:
		tem = inl.split()
		ch = tem[1].split(':')[0]
		region_left = int(tem[1].split(':')[1].split('-')[0])-1   ### region left
		region_right = int(tem[1].split(':')[1].split('-')[1])   ### region right
		if region_right > len(S[ch]):
			region_right = len(S[ch])
		E_left = max(region_left - l,0)  ### extend left
		E_right = min(region_right + l,len(S[ch]))   ### extend right
		seq_left = S[ch][E_left:region_left]  ### left sequence
		seq_right = S[ch][region_right:E_right]   ### right sequence
		GC = GC_content(seq_left + seq_right)
		gene_length = 0
		repeat_length = 0
		tandem_length = 0
		E_length = LWN(seq_left) + LWN(seq_right)
		if E_length != 0:
			### get total gene length for neighbor sequences
			if ch in G:
				for gene in G[ch].keys():
					if (region_left-G[ch][gene][0])*(E_left-G[ch][gene][1]) < 0 or E_left-G[ch][gene][0]==0 or region_left-G[ch][gene][1]==0:
						length = LWN(S[ch][G[ch][gene][0]:G[ch][gene][1]])
						distance = min(LWN(S[ch][min(region_left,G[ch][gene][0]):max(region_left,G[ch][gene][0])]),LWN(S[ch][min(E_left,G[ch][gene][1]):max(E_left,G[ch][gene][1])]),length,LWN(seq_left))
						gene_length += distance
						if gene in T:
							tandem_length += distance
					if (E_right-G[ch][gene][0])*(region_right-G[ch][gene][1]) < 0 or region_right-G[ch][gene][0]==0 or E_right-G[ch][gene][1]==0:
						length = LWN(S[ch][G[ch][gene][0]:G[ch][gene][1]])
						distance = min(LWN(S[ch][min(G[ch][gene][0],E_right):max(G[ch][gene][0],E_right)]),LWN(S[ch][min(G[ch][gene][1],region_right):max(G[ch][gene][1],region_right)]),length,LWN(seq_right))
						gene_length += distance
						if gene in T:
							tandem_length += distance
			### get total repeat length for neighbor sequences
			if ch in R:
				for gene in R[ch].keys():
					if (region_left-R[ch][gene][0])*(E_left-R[ch][gene][1]) < 0 or E_left-R[ch][gene][0]==0 or region_left-R[ch][gene][1]==0:
						length = LWN(S[ch][R[ch][gene][0]:R[ch][gene][1]])
						distance = min(LWN(S[ch][min(region_left,R[ch][gene][0]):max(region_left,R[ch][gene][0])]),LWN(S[ch][min(E_left,R[ch][gene][1]):max(E_left,R[ch][gene][1])]),length,LWN(seq_left))
						repeat_length += distance
					if (E_right-R[ch][gene][0])*(region_right-R[ch][gene][1]) < 0 or region_right-R[ch][gene][0]==0 or E_right-R[ch][gene][1]==0:
						length = LWN(S[ch][R[ch][gene][0]:R[ch][gene][1]])
						distance = min(LWN(S[ch][min(R[ch][gene][0],E_right):max(R[ch][gene][0],E_right)]),LWN(S[ch][min(R[ch][gene][1],region_right):max(R[ch][gene][1],region_right)]),length,LWN(seq_right))
						repeat_length += distance
			out.write('%s\t%s:%s_%s\t%s\t%s\t%s\t%s\t%s\n'%(type,ch,region_left,region_right,l,GC,float(gene_length)/E_length,float(repeat_length)/E_length,float(tandem_length)/E_length))
			out.flush()

out.close()


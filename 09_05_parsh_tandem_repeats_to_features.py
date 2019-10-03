import sys,os

genome_file = sys.argv[1] #Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fna.longest.mod.fa
CNV_file = sys.argv[2] 
SSR_file = sys.argv[3] #Tandem_repeats_whole_genome.txt

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

tr = open(SSR_file,'r').readlines()[1:]
### Note that some Tandem_repeats would be overlapping with each others
T_pre = {} ### for Tandem_repeats
for inl in tr:
	tem = inl.split('\t')
	chr = tem[0]
	left = int(tem[1])
	right = int(tem[2])
	gene = '%s_%s_%s'%(chr,left,right)
	if chr not in T_pre:
		T_pre[chr] = {}
	T_pre[chr][left] = [gene,right]

nn = 0
total_nn = 0
for ch in sorted(T_pre.keys()):
	for left in sorted(T_pre[ch].keys()):
		total_nn += len(T_pre[ch][left])

while total_nn != nn:
	nn = total_nn
	for ch in sorted(T_pre.keys()):
		left = sorted(T_pre[ch].keys())
		x = 0
		while x < len(left)-1:
			y = x + 1
			while y < len(left) and T_pre[ch][left[x]][1] >= left[y]:
				T_pre[ch][left[x]][1] = T_pre[ch][left[y]][1]
				T_pre[ch].pop(left[y])
				left.pop(y)
				y += 1
			x = y
	total_nn = 0
	for ch in sorted(T_pre.keys()):
		for left in sorted(T_pre[ch].keys()):
			total_nn += len(T_pre[ch][left])

G = {} ### for tandem repeat
for ch in sorted(T_pre.keys()):
	if ch not in G:
		G[ch] = {}
	for left in sorted(T_pre[ch].keys()):
		G[ch][T_pre[ch][left][0]] = [left,T_pre[ch][left][1]]


inp = open(CNV_file,'r').readlines()
out = open(CNV_file + '_Tandem_repeats_density','w')
out.write('Region\tRD\tTandem_repeats_density\n')
for inl in inp:
	tem = inl.strip().split('\t')
	ch = tem[0].split(':')[0]
	region_left = int(tem[0].split(':')[1].split('-')[0])
	region_right = int(tem[0].split(':')[1].split('-')[1])
	gene_length = 0
	region_length = LWN(S[ch][region_left:region_right])
	if region_length == 0:
		out.write(inl.strip() + '\t%s\n'%(0))
	else:
		if ch in G:
			for gene in G[ch].keys():
				if (region_left-G[ch][gene][1])*(region_right-G[ch][gene][0]) < 0 or region_left-G[ch][gene][0]==0 or region_right-G[ch][gene][1]==0:
					length = LWN(S[ch][G[ch][gene][0]:G[ch][gene][1]])
					distance = min(LWN(S[ch][min(region_left,G[ch][gene][1]):max(region_left,G[ch][gene][1])]),LWN(S[ch][min(region_right,G[ch][gene][0]):max(region_right,G[ch][gene][0])]),length,region_length)
					gene_length += distance
		out.write(inl.strip() + '\t%s\n'%(float(gene_length)/region_length))
		out.flush()

out.close()

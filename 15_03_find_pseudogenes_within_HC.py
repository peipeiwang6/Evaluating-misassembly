import sys,os
HC_file = sys.argv[1] #observed or simulated HC regions
gff_file = sys.argv[2] #sl_ps.gff
genome_file = sys.argv[3] #Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fna.longest.mod.fa
output = sys.argv[4]

cnv = open(HC_file,'r').readlines()
gff = open(gff_file,'r').readlines()
genome = open(genome_file,'r').readlines()
out = open(output,'w')

GENOME = {}
x = 0
while x < len(genome)/2:
	chr = genome[2*x].strip()[5:]
	seq = genome[2*x+1]
	GENOME[chr] = seq
	x += 1

def region_length(seq):
	n = 0
	for s in seq:
		if s != 'N' and s != 'n':
			n += 1
	return n

G = {}
x = 0
while x < len(gff):
	gene = gff[x].split('\t')[1]
	ch = gff[x].split('\t')[0]
	if ch != 'NC_007898.3':
		if ch not in G:
			G[ch] = {}
		G[ch][gene] = [int(gff[x].split('\t')[2]),int(gff[x].split('\t')[3].strip()),'pseudo']
	x += 1

out.write('CNVs\tGene\ttype\tChr\tGene_left\tGene_right\tProportion_of_gene\tRegion_length_without_Ns\tRegion_length\tGene_length\n')
Result = {}
for inl in cnv:
	if sys.argv[2] in inl:
		ch = inl.split('\t')[1].split(':')[0]
		if ch != 'NC_007898.3':
			if ch in G:
				left = int(inl.split('\t')[1].split(':')[1].split('-')[0])
				right = int(inl.split('\t')[1].split(':')[1].split('-')[1])
				for gene in G[ch].keys():
					if (right-G[ch][gene][0])*(left-G[ch][gene][1]) < 0 or left-G[ch][gene][0]==0 or right-G[ch][gene][1]==0:
						length = G[ch][gene][1] - G[ch][gene][0] + 1
						distance = min(abs(right-G[ch][gene][0]),abs(left-G[ch][gene][1]),length, right-left+1)
						if float(distance)/float(length) > 1:
							overlap = 1
						else:
							overlap = str(float(distance)/float(length))
						out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(inl.split('\t')[1],gene,G[ch][gene][2],ch,G[ch][gene][0],G[ch][gene][1],overlap,region_length(GENOME[ch][left-1:right]),right-left+1,length))
						Result[inl.split('\t')[1]] = 1
				if inl.split('\t')[1] not in Result:
					out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(inl.split('\t')[1],'No_gene','No_gene','No_gene',0,0,0,0,right-left+1,0))
			else:
				out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(inl.split('\t')[1],'No_gene','No_gene','No_gene',0,0,0,0,right-left+1,0))

out.close()
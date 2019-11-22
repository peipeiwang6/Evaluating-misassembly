import sys,os
## SL4.0
inp = open('S_lycopersicum_chromosomes.4.00.fa','r')
out = open('S_lycopersicum_chromosomes.4.00_mod.fa','w')
inl = inp.readline()
while inl:
	if inl.startswith('>'):
		out.write('\n'+inl)
	else:
		out.write(inl.strip())
	inl = inp.readline()
	
out.close()


inp = open('S_lycopersicum_chromosomes.4.00_mod.fa','r')
out = open('S_lycopersicum_chromosomes.4.00_merged.fa','w')
out2 = open('Sly_v4_merged_chr_length.fa','w')
out.write('>SL4.0\n')
SEQ = {}
L = {}
inl = inp.readline()
while inl:
	if inl.startswith('>'):
		chr = inl.strip()[1:]
		L[chr] = 1
		SEQ[chr] = ''
	else:
		L[chr] = len(inl.strip())
		SEQ[chr] = inl.strip()
	inl = inp.readline()

for chrom in ['SL4.0ch01','SL4.0ch02','SL4.0ch03','SL4.0ch04','SL4.0ch05','SL4.0ch06','SL4.0ch07','SL4.0ch08','SL4.0ch09','SL4.0ch10','SL4.0ch11','SL4.0ch12','SL4.0ch00']:
	out.write(SEQ[chrom])
	out2.write('%s\t%s\n'%(chrom,len(SEQ[chrom])))
	
out.close()
out2.close()

# SL2.5
inp = open('Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fna.longest.mod.fa','r')
out = open('SL2.5_merged.fa','w')
out2 = open('SL2.5_merged_chr_length.fa','w')
out.write('>SL2.5\n')
SEQ = {}
L = {}
inl = inp.readline()
while inl:
	if inl.startswith('>'):
		chr = inl.strip()[1:]
		L[chr] = 1
		SEQ[chr] = ''
	else:
		L[chr] = len(inl.strip())
		SEQ[chr] = inl.strip()
	inl = inp.readline()

chr = sorted(L.keys())
chr_ordered = chr[1:13]
chr_ordered.append(chr[0])
for j in range(13,len(chr)):
	chr_ordered.append(chr[j])

for chrom in chr_ordered:
	out.write(SEQ[chrom])
	out2.write('%s\t%s\n'%(chrom,len(SEQ[chrom])))
	
out.close()
out2.close()
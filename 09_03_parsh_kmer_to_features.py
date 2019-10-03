import sys,os
import itertools
from itertools import product

genome_file = sys.argv[1] #Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fna.longest.mod.fa
CNV_file = sys.argv[2] 
kmer_file = sys.argv[3] 
k = int(sys.argv[4]) 

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


nucleotide = 'ATCG'
if k == 1:
	N = ['A','T','C','G']
if k > 1:
	N = []
	nn = itertools.product('ATCG', repeat=k)
	for n in nn:
		N.append(''.join(n))

N = sorted(N)

kmer = open(kmer_file,'r').readlines()
D = {}
for inl in kmer:
	tem = inl.strip().split('\t')
	region = tem[0]
	nn = tem[2]
	num = float(tem[3])
	if nn not in D:
		D[nn] = {}
	D[nn][region] = num


inp = open(CNV_file,'r').readlines()
out = open(CNV_file + '_Kmer_%s'%k,'w')
title = 'Region\tRD'
for n in N:
	title = title + '\t' + n

out.write(title + '\n')

for inl in inp:
	tem = inl.strip().split('\t')
	ch = tem[0].split(':')[0]
	region = tem[0]
	left = int(tem[0].split(':')[1].split('-')[0])
	right = int(tem[0].split(':')[1].split('-')[1])
	length = LWN(S[ch][left:right])
	res = inl.strip()
	for nn in N:
		if nn in D:
			if region in D[nn]:
				if length == 0:
					res = res + '\t0'
				else:
					res = res + '\t%s'%(D[nn][region]/length)
			else:
				res = res + '\t0'
		else:
			res = res + '\t0'
	out.write(res+'\n')

out.close()
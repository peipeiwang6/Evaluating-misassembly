'''
input 1: path to save your output files
input 2: HC regions with high confidence
input 3: LC regions with high confidence
input 4: length of chromosomes
'''
import sys,os
import random
import decimal

path = sys.argv[1]
HC_file = sys.argv[2]
LC_file = sys.argv[3]
Length_chr = sys.argv[4]

duplication = open(HC_file,'r').readlines()
for inl in duplication:
	region = inl.split('\t')[1]
	Length = int(decimal.Decimal(inl.split('\t')[2]))
	D['HC__'+region] = Length

deletion = open(LC_file,'r').readlines()
D = {}
for inl in deletion:
	region = inl.split('\t')[1]
	Length = int(decimal.Decimal(inl.split('\t')[2]))
	D['LC__'+region] = Length

key_sorted = sorted(D, key=D.get,reverse=True)  ### sort the HC and LC regions, for the longest to the shortest

chromo = open(Length_chr,'r').readlines()
C = {}
L = {}
for inl in chromo:
	l = int(inl.strip().split('\t')[1])
	chr = inl.split('\t')[0]
	C[chr] = l

background = open(Length_chr,'r').readlines()
L = {}
for inl in background:
	b = int(inl.strip().split('\t')[1])
	chr = inl.split('\t')[0]
	loc = range(1,b,100)
	for ll in loc:
		L['%s-%s'%(chr,ll)] = ll

def LOCUS(locus,x):
	for aa in x:
		L.pop('%s-%s'%(locus.split('-')[0],aa))

def LOCUS_in_L(locus,x):  ### whether all the locus still in the regions
	for aa in x:
		if '%s-%s'%(locus.split('-')[0],aa) in L:
			res = 'T'
		else:
			res = 'F'
			break
	return res
	
for i in range(1,1001):
	out = open(path + 'Randomly_choose_HC_LC_result_%s.txt'%(i),'w')
	for r in key_sorted:
		n = 0
		while n == 0:
			y = 0
			locus = random.choice(list(L.keys()))
			length = int(locus.split('-')[1]) + D[r] -1
			if length <= C[locus.split('-')[0]] and LOCUS_in_L(locus,range(L[locus],length,100))=='T':
				out1 = r.split('__')[0].strip()
				out2 = locus.split('-')[0].strip()
				out3 = locus.split('-')[1].strip()
				out4 = length
				out5 = D[r]
				out.write('%s\t%s:%s-%s\t%s\n'%(out1,out2,out3,out4,out5))
				print(out1+"\t"+out2+':'+out3+'-'+str(out4)+'\t'+str(out5))
				LOCUS(locus,range(L[locus],length,100))
				n = 1
			else:
				n = 0
	out.close()

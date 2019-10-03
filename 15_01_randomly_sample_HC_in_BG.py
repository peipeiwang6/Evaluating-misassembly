import sys,os
import random
import decimal

HC_file = sys.argv[1] # 20180225_Duplication_region_with_0.7_1.76_coverage.txt
BG_file = sys.argv[2] # 20180225_Background_region_with_0.7_1.76_coverage.txt_modified
chromosome_length = sys.argv[3] #Sly_chromosome_length.txt
path = sys.argv[4] # path to save random selected HC regions

dup = open(HC_file,'r').readlines()
D = {}
for inl in dup:
	if 'duplication' in inl:
		region = inl.split('\t')[1]
		length = int(decimal.Decimal(inl.split('\t')[2]))
		D[region] = length

key_sorted = sorted(D, key=D.get,reverse=True)  ### sort the duplication region, for the longest to the shortest

chromo = open(chromosome_length,'r').readlines()
C = {}
L = {}
for inl in chromo:
	l = int(inl.strip().split('\t')[1])
	chr = inl.split('\t')[0]
	C[chr] = l

background = open(BG_file,'r').readlines()
L = {}
for inl in background:
	l = int(inl.split('\t')[1].split(':')[1].split('-')[0])
	r = int(inl.split('\t')[1].split(':')[1].split('-')[1])
	chr = inl.split('\t')[1].split(':')[0]
	loc = range(l,r,100)
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

for i in range(1,10001,1):
	out = open(path + 'Randomly_HC_region_%s.txt'%(i),'w')
	for r in key_sorted:
		n = 0
		while n == 0:
			y = 0
			locus = random.choice(L.keys())
			length = int(locus.split('-')[1]) + D[r]
			if length <= C[locus.split('-')[0]] and LOCUS_in_L(locus,range(L[locus],length+100,100))=='T':
				out.write('duplication\t%s:%s-%s\t%s\n'%(locus.split('-')[0],locus.split('-')[1],length,D[r]))
				LOCUS(locus,range(L[locus],length+100,100))
				print locus
				n = 1
			else:
				n = 0
	out.close()


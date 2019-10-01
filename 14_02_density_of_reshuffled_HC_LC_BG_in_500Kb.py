'''
input 1: files with length of chromosomes
input 2: correlation of chromosome names between genome assembly version NCBI v2.5 and Solanaceae Genomics Network ITAG2.4 release
input 3: path with the reshuffled HC/LC regions
'''
import sys,os
Length_chr = sys.argv[1]
Dic = sys.argv[2]
path = sys.argv[3]

chr = open(Length_chr,'r').readlines()
CH = {}
for inl in chr:
	if inl.startswith('NC') and not inl.startswith('NC_007898.3'):
		CH[inl.split('\t')[0]] = int(inl.split('\t')[1].strip())
	
dic = open(Dic,'r').readlines()
D = {}
for inl in dic:
	D[inl.split('\t')[0]] = inl.split('\t')[1].strip()

for n in range(1,1001):	
	file = path + 'Randomly_choose_HC_LC_result_%s.txt'%n
	inp = open(file,'r').readlines()
	### get background regions
	L = {}
	H = {}
	B = {}
	All = {}
	for inl in inp:
		try:
			chromosome = inl.split('\t')[1].split(':')[0]
			if chromosome.startswith('NC') and not chromosome.startswith('NC_007898.3'):
				region_left = int(inl.split('\t')[1].split(':')[1].split('-')[0])
				region_right = int(inl.split('\t')[1].split(':')[1].split('-')[1])
				type = inl.split('\t')[0]
				if chromosome not in All:
					All[chromosome] = {}
				All[chromosome][region_left] = region_right
				if type == 'LC':
					if chromosome not in L:
						L[chromosome] = {}
					L[chromosome][region_left] = region_right
				if type == 'HC':
					if chromosome not in H:
						H[chromosome] = {}
					H[chromosome][region_left] = region_right
		except:
			print('Error with %s'%inl)
			
	for chromosome in All:
		if chromosome not in B:
			B[chromosome] = {}
		sorted_left = sorted(All[chromosome].keys())
		if sorted_left[0] != 1:
			B[chromosome][1] = sorted_left[0]-1
		for i in range(1,len(sorted_left)):
			B[chromosome][All[chromosome][sorted_left[i-1]] + 1] = sorted_left[i] -1
		B[chromosome][All[chromosome][sorted_left[-1]] + 1] = CH[chromosome]
		
		
	out = open(file + '_density','w')
	out.write('Chr\tNumber_of_bin\tHC_density\tbackground_density\tLC_density\n')
	R = {}
	for ch in CH.keys():
		if ch not in R:
			R[ch] = {}
		for i in range(1,CH[ch]/500000+1):
			left = (i-1)*500000+1
			right = i *500000
			R[ch][i] = {}
			R[ch][i]['HC'] = 0
			R[ch][i]['LC'] = 0
			R[ch][i]['background'] = 0
			for region_left in H[ch]:
				if (right - region_left)*(left - H[ch][region_left]) <= 0:
					overlapping = min(abs(H[ch][region_left] - left + 1), abs(right-region_left+1),abs(right-left+1),abs(H[ch][region_left] - region_left + 1))
					R[ch][i]['HC'] += overlapping
			for region_left in L[ch]:
				if (right - region_left)*(left - L[ch][region_left]) <= 0:
					overlapping = min(abs(L[ch][region_left] - left + 1), abs(right-region_left+1),abs(right-left+1),abs(L[ch][region_left] - region_left + 1))
					R[ch][i]['LC'] += overlapping
			for region_left in B[ch]:
				if (right - region_left)*(left - B[ch][region_left]) <= 0:
					overlapping = min(abs(B[ch][region_left] - left + 1), abs(right-region_left+1),abs(right-left+1),abs(B[ch][region_left] - region_left + 1))
					R[ch][i]['background'] += overlapping
					
	for ch in R:
		for i in R[ch]:
			out.write('%s\t%s\t%s\t%s\t%s\n'%(ch,i,float(R[ch][i]['HC'])/500000,float(R[ch][i]['LC'])/500000,float(R[ch][i]['background'])/500000))
			
	out.close()
	
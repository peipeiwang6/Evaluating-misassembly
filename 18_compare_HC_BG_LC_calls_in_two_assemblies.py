import sys,os
chr_old = open('SL2.5_merged_chr_length.txt','r').readlines()
chr_new = open('Sly_v4_merged_chr_length.txt','r').readlines()
L_old = {}
L_new = {}
p_old = 0
p_new = 0
for inl in chr_old:
	L_old[inl.split('\t')[0]] = [p_old+1,p_old + int(inl.split('\t')[1].strip())]
	p_old = p_old + int(inl.split('\t')[1].strip())
	
for inl in chr_new:
	L_new[inl.split('\t')[0]] = [p_new+1,p_new + int(inl.split('\t')[1].strip())]
	p_new = p_new + int(inl.split('\t')[1].strip())


def overlapping_perc(left1,right1,left2,right2):
	if (left1-right2)*(right1-left2) < 0:
		overlap = min(abs(left1-right2)+1, abs(right1-left2)+1,abs(left1-right1) + 1,abs(left2-right2) + 1)
		return(overlap/(abs(left1-right1) + 1),overlap/(abs(left2-right2) + 1),overlap)
	else:
		return(0,0,0)
	

correponding = open('All_chr_corresponding_between_two_assemblies.txt','r').readlines()
out = open('All_chr_corresponding_between_two_assemblies_true_chr_pos.txt','w')
for inl in correponding:
	tem = inl.split('\t')
	left1 = min(int(tem[1]),int(tem[2]))
	right1 = max(int(tem[1]),int(tem[2]))
	OL = {}
	for chr in L_old:
		left2 = L_old[chr][0]
		right2 = L_old[chr][1]
		if overlapping(left1,right1,left2,right2) > 0:
			OL[overlapping(left1,right1,left2,right2)] = chr
	chr = OL[max(OL.keys())]
	tem[0] = chr
	tem[1] = '%s'%(int(tem[1]) - L_old[chr][0] + 1)
	tem[2] = '%s'%(int(tem[2]) - L_old[chr][0] + 1)
	left1 = min(int(tem[4]),int(tem[5]))
	right1 = max(int(tem[4]),int(tem[5]))
	OL = {}
	for chr in L_new:
		left2 = L_new[chr][0]
		right2 = L_new[chr][1]
		if overlapping(left1,right1,left2,right2) > 0:
			OL[overlapping(left1,right1,left2,right2)] = chr
	chr = OL[max(OL.keys())]
	tem[3] = chr
	tem[4] = '%s'%(int(tem[4]) - L_new[chr][0] + 1)
	tem[5] = '%s'%(int(tem[5]) - L_new[chr][0] + 1)
	out.write('\t'.join(tem))

out.close()	


HC_old = open('20180225_Duplication_region_with_0.7_1.76_coverage.txt','r').readlines()
HC_new = open('HC_new_assembly.txt','r').readlines()

Old_HC = {}
New_HC = {}
for inl in HC_old:
	tem = inl.split('\t')
	chr = tem[1].split(':')[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	if chr not in Old_HC:
		Old_HC[chr]=[]
	Old_HC[chr].append([left,right,tem[2],tem[3].strip()])
	
for inl in HC_new:
	tem = inl.split('\t')
	chr = tem[1].split(':')[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	if chr not in New_HC:
		New_HC[chr]=[]
	New_HC[chr].append([left,right,tem[2],tem[3].strip()])
	
correponding = open('All_chr_corresponding_between_two_assemblies_true_chr_pos.txt','r').readlines()
Corr = []
for inl in correponding:
	if 'reverse' not in inl:
		tem = inl.split('\t')
		Corr.append([tem[0],int(tem[1]),int(tem[2]),tem[3],int(tem[4]),int(tem[5])])

out = open('Corresponding_HC_BC_regions.txt','w')	
for corr in Corr:
	chr_old = corr[0]
	chr_new = corr[3]
	result = '%s\t%s\t%s\t%s\t%s\t%s'%(corr[0],corr[1],corr[2],corr[3],corr[4],corr[5])
	if chr_old in Old_HC:
		for loc in Old_HC[chr_old]:
			res = overlapping_perc(min(corr[1],corr[2]),max(corr[1],corr[2]),loc[0],loc[1])
			if res[2] != 0:
				out.write(result + '\tOld\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tHC'%(loc[0],loc[1],loc[2],loc[3],res[0],res[1],res[2]) + '\n')
	if chr_new in New_HC:
		for loc in New_HC[chr_new]:
			res = overlapping_perc(min(corr[4],corr[5]),max(corr[4],corr[5]),loc[0],loc[1])
			if res[2] != 0:
				out.write(result + '\tNew\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tHC'%(loc[0],loc[1],loc[2],loc[3],res[0],res[1],res[2]) + '\n')
	

BC_old = open('20180225_Background_region_with_0.7_1.76_coverage.txt_modified','r').readlines()
BC_new = open('BC_new_assembly_filted.txt','r').readlines()

Old_BC = {}
New_BC = {}
for inl in BC_old:
	tem = inl.split('\t')
	chr = tem[1].split(':')[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	if chr not in Old_BC:
		Old_BC[chr]=[]
	Old_BC[chr].append([left,right,tem[2],tem[3].strip()])
	
for inl in BC_new:
	tem = inl.split('\t')
	chr = tem[1].split(':')[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	if chr not in New_BC:
		New_BC[chr]=[]
	New_BC[chr].append([left,right,tem[2],tem[3].strip()])
	
for corr in Corr:
	chr_old = corr[0]
	chr_new = corr[3]
	result = '%s\t%s\t%s\t%s\t%s\t%s'%(corr[0],corr[1],corr[2],corr[3],corr[4],corr[5])
	if chr_old in Old_BC:
		for loc in Old_BC[chr_old]:
			res = overlapping_perc(min(corr[1],corr[2]),max(corr[1],corr[2]),loc[0],loc[1])
			if res[2] != 0:
				out.write(result + '\tOld\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tBC'%(loc[0],loc[1],loc[2],loc[3],res[0],res[1],res[2]) + '\n')
	if chr_new in New_BC:
		for loc in New_BC[chr_new]:
			res = overlapping_perc(min(corr[4],corr[5]),max(corr[4],corr[5]),loc[0],loc[1])
			if res[2] != 0:
				out.write(result + '\tNew\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tBC'%(loc[0],loc[1],loc[2],loc[3],res[0],res[1],res[2]) + '\n')
	
LC_old = open('20180225_Deletion_region_with_0.7_1.76_coverage.txt','r').readlines()
LC_new = open('LC_new_assembly.txt','r').readlines()
Old_LC = {}
New_LC = {}
for inl in LC_old:
	tem = inl.split('\t')
	chr = tem[1].split(':')[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	if chr not in Old_LC:
		Old_LC[chr]=[]
	Old_LC[chr].append([left,right,tem[2],tem[3].strip()])
	
for inl in LC_new:
	tem = inl.split('\t')
	chr = tem[1].split(':')[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	if chr not in New_LC:
		New_LC[chr]=[]
	New_LC[chr].append([left,right,tem[2],tem[3].strip()])
	
for corr in Corr:
	chr_old = corr[0]
	chr_new = corr[3]
	result = '%s\t%s\t%s\t%s\t%s\t%s'%(corr[0],corr[1],corr[2],corr[3],corr[4],corr[5])
	if chr_old in Old_LC:
		for loc in Old_LC[chr_old]:
			res = overlapping_perc(min(corr[1],corr[2]),max(corr[1],corr[2]),loc[0],loc[1])
			if res[2] != 0:
				out.write(result + '\tOld\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tLC'%(loc[0],loc[1],loc[2],loc[3],res[0],res[1],res[2]) + '\n')
	if chr_new in New_LC:
		for loc in New_LC[chr_new]:
			res = overlapping_perc(min(corr[4],corr[5]),max(corr[4],corr[5]),loc[0],loc[1])
			if res[2] != 0:
				out.write(result + '\tNew\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tLC'%(loc[0],loc[1],loc[2],loc[3],res[0],res[1],res[2]) + '\n')
	
	
out.close()


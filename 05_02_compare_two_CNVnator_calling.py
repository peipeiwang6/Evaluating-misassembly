from __future__ import division 
import os,sys
import random
import time
import numpy as np
fs=open(sys.argv[1],"r") ## the CNVnator calling results using simulated reads
fo=open(sys.argv[2],'r') ## the original, or analysis CNVnator calling results, using reads of dataset2
fc=open(sys.argv[3],"w") ## output file
fo_d={}
for oline in fo.readlines():
	tem=oline.strip().split("\t")
	region = tem[1]
	ty = tem[0]
	ch = tem[1].split(":")[0]
	left = int(tem[1].split(':')[1].split('-')[0])
	right = int(tem[1].split(':')[1].split('-')[1])
	fo_d[region] = [ch,left,right,ty]
fs_d = {}
for sline in fs.readlines():
	tem = sline.strip().split("\t")
        region = tem[1]
	ty = tem[0]
        ch = tem[1].split(":")[0]
        left = int(tem[1].split(':')[1].split('-')[0])
        right = int(tem[1].split(':')[1].split('-')[1])
        fs_d[region] = [ch,left,right,ty]
def overlap(ch,left,right,Fsd,ty):
	res = []
	for region in Fsd.keys():
		if Fsd[region][0] == ch and ty == Fsd[region][3]:
			if (Fsd[region][1] >= left and Fsd[region][1] <= right) or (Fsd[region][2] >= left and Fsd[region][2] <= right) or (Fsd[region][1] < left and Fsd[region][2] > right):
				res.append(region)
	return res
Results = []
for region in fo_d.keys():
	target = overlap(fo_d[region][0],fo_d[region][1],fo_d[region][2],fs_d,fo_d[region][3])
	left_r=int(region.split(":")[1].split("-")[0])
	right_r=int(region.split(":")[1].split("-")[1])
	tmplen=0
	for i in range(len(target)):
		left_t=int(target[i].split(":")[1].split("-")[0])
		right_t=int(target[i].split(":")[1].split("-")[1])
		if left_t>=left_r and left_t<=right_r and right_t>right_r:
			length_TP=right_r-left_t+1
		elif left_t<left_r and right_t>=left_r and right_t <= right_r:
			length_TP=right_t-left_r+1
		elif left_t>=left_r and right_t<=right_r:
			length_TP=right_t-left_t +1
		elif left_t <= left_r and right_t >=right_r and ((right_t-left_t)!=(right_r-left_r)):
			length_TP=right_r-left_r+1
		tmplen=tmplen+length_TP
	sensitivity= int(tmplen) / int((right_r-left_r+1))
	Results.append([fo_d[region][3],region,int(tmplen),right_r-left_r+1,str(sensitivity)])
	
Results = np.array(Results)
total_num = len(Results[:,0])
total_dup = len(Results[Results[:,0]=="duplication"])
total_del = len(Results[Results[:,0]=="deletion"])
total_det_array = Results[Results[:,4]!=str(0.0)]
total_det = len(total_det_array)
total_det_dup = len(total_det_array[total_det_array[:,0]=="duplication"])
total_det_del = len(total_det_array[total_det_array[:,0]=="deletion"])
total_no_det_array = Results[Results[:,4]==str(0.0)]
total_no_det = len(total_no_det_array)
total_no_det_dup = len(total_no_det_array[total_no_det_array[:,0]=="duplication"])
total_no_det_del = len(total_no_det_array[total_no_det_array[:,0]=="deletion"])
TP_nucleotide_array = Results[:,2].astype(int)
Total_nucleotide_array = Results[:,3].astype(int)
TP_nucleotide = np.sum(TP_nucleotide_array)
Total_nucleotide = np.sum(Results[:,3].astype(int))
Nuc_sensitivity = TP_nucleotide/Total_nucleotide
Fregment_sensitivity = total_det/total_num
fc.write("Nuc_Precision\tFregment_Precision\n")
fc.write("%s\t%s\n"%(Nuc_sensitivity,Fregment_sensitivity))

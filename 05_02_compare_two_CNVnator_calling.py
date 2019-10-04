from __future__ import division 
import os,sys
import random
import time
import numpy as np
fs_cnv=open(sys.argv[1],"r") ## the CNVnator calling results using simulated reads
fo_cnv=open(sys.argv[2],'r') ## the original, or analysis CNVnator calling results, using reads of dataset2
fc=open(sys.argv[3],"w") ## output file
def overlap(ch,left,right,Fsd,ty):
	res = []
	for region in Fsd.keys():
		if Fsd[region][0] == ch and ty == Fsd[region][3]:
			if (Fsd[region][1] >= left and Fsd[region][1] <= right) or (Fsd[region][2] >= left and Fsd[region][2] <= right) or (Fsd[region][1] < left and Fsd[region][2] > right):
				res.append(region)
	return res
def readcnv(cnv):
	fd={}
	for oline in cnv.readlines():
		temp=oline.strip().split("\t")
		region = temp[1]
		ty = temp[0]
		ch = temp[1].split(":")[0]
		left = int(temp[1].split(':')[1].split('-')[0])
		right = int(temp[1].split(':')[1].split('-')[1])
		fd[region] = [ch,left,right,ty]
	return(fd)
fod = readcnv(fo_cnv)
fsd = readcnv(fs_cnv)
def compare(fo,fs):
	Results = []
	for region in fo.keys():
		target = overlap(fo[region][0],fo[region][1],fo[region][2],fs,fo[region][3])
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
		Results.append([fo[region][3],region,int(tmplen),right_r-left_r+1,str(sensitivity)])
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
	return(float(Nuc_sensitivity),float(Fregment_sensitivity))
Recall_nuc,Recall_freg=tuple(compare(fsd,fod))
Precision_nuc,Precision_freg=tuple(compare(fod,fsd))
F1_nuc = 2*(Recall_nuc*Precision_nuc)/(Recall_nuc+Precision_nuc)
F1_freg = 2*(Recall_freg*Precision_freg)/(Recall_freg+Precision_freg)
fc.write("Nuc_F1\tFregment_F1\n")
fc.write("%s\t%s\n"%(F1_nuc,F1_freg))

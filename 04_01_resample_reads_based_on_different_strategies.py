import os,sys
import argparse
import random
import time
from Bio import SeqIO

CNV_file = sys.argv[1] # the copy number file
genome = sys.argv[2] # genome file
strategies = sys.argv[3] # strategies to resample reads

fasta_sequences = SeqIO.parse(open(genome),'fasta')
# separate chromosomes into individual files, endswith .fa 
for seqs in fasta_sequences:
	tem=open(seqs.id+".fa", "w")
	SeqIO.write([seqs], tem, "fasta")
	
fcnv=open(CNV_file)

if strategies=='012':
	fsimu=open("Sly_simulation_Discretized_to_0_1_2_"+str(time.strftime("%Y-%m-%d_%H-%M-%S",time.localtime()))+".fq","a")
	for line in fcnv.readlines():
		cline=line.strip().split("\t")
		if cline[0]=="deletion":
			readsnum=0
		elif cline[0] == 'duplication':
			readsnum=2
		else:
			readsnum=1
		cline.append(round(int(float(cline[2]))*46.3*readsnum/101))  ### here, 46.3 is the read coverage of datasets, 101 is the read length
		cline[1]=cline[1].split(":")
		cline[1][1]=cline[1][1].split("-")
		cline.append([])
		if cline[-2] >0:
			chromosomes=open(str(cline[1][0]).strip()+".fa")
			tem=""
			for gen in chromosomes.readlines():
					if gen.startswith(">"):
							tem=""
					else:
							tem=tem+gen.strip()
			for n in range(1,int(float((cline[-2])))):
				rd=random.randint(int(float(cline[1][1][0])),int(float(cline[1][1][1])))
				temreads=tem[(rd-51):(rd+50)].strip()
				if "N" not in temreads and temreads != "":
					fsimu.write("@Simu_"+str(cline[1][0]).strip()+"_"+str(n)+"_"+str(rd)+"\n"+temreads+"\n"+"+"+"\n"+"h"*len(temreads)+"\n")

if strategies=='rounded':
	fsimu=open("Sly_simulation_rounded_"+str(time.strftime("%Y-%m-%d_%H-%M-%S",time.localtime()))+".fq","a")
	for line in fcnv.readlines():
		cline=line.strip().split("\t")
		cline.append(int(round((int(float(cline[2]))*46.3*round(float(cline[3]))/101))))
		cline[1]=cline[1].split(":")
		cline[1][1]=cline[1][1].split("-")
		cline.append([])
		if cline[-2] >0:
			chromosomes=open(str(cline[1][0]).strip()+".fa")
			tem=""
			for gen in chromosomes.readlines():
					if gen.startswith(">"):
							tem=""
					else:
							tem=tem+gen.strip()
			for n in range(1,int(float((cline[-2])))):
				rd=random.randint(int(float(cline[1][1][0])),int(float(cline[1][1][1])))
				temreads=tem[(rd-51):(rd+50)].strip()
				if "N" not in temreads and temreads != "":
					fsimu.write("@Simu_"+str(cline[1][0]).strip()+"_"+str(n)+"_"+str(rd)+"\n"+temreads+"\n"+"+"+"\n"+"h"*len(temreads)+"\n")

if strategies=='analysis':
	fsimu=open("Sly_simulation_analysis_"+str(time.strftime("%Y-%m-%d_%H-%M-%S",time.localtime()))+".fq","a")
	for line in fcnv.readlines():
		cline=line.strip().split("\t")
		cline.append(int(round((int(float(cline[2]))*46.3*float(cline[3])/101))))
		cline[1]=cline[1].split(":")
		cline[1][1]=cline[1][1].split("-")
		cline.append([])
		if cline[-2] >0:
			chromosomes=open(str(cline[1][0]).strip()+".fa")
			tem=""
			for gen in chromosomes.readlines():
					if gen.startswith(">"):
							tem=""
					else:
							tem=tem+gen.strip()
			for n in range(1,int(float((cline[-2])))):
				rd=random.randint(int(float(cline[1][1][0])),int(float(cline[1][1][1])))
				temreads=tem[(rd-51):(rd+50)].strip()
				if "N" not in temreads and temreads != "":
					fsimu.write("@Simu_"+str(cline[1][0]).strip()+"_"+str(n)+"_"+str(rd)+"\n"+temreads+"\n"+"+"+"\n"+"h"*len(temreads)+"\n")

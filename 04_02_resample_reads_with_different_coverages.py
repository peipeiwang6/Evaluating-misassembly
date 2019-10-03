import sys,os
import random

number = int(sys.argv[1]) ### number of read coverage you want to resample 
output1 = open('Resample_%s_1.fastq'%number,'w')
output2 = open('Resample_%s_2.fastq'%number,'w')
D = []
genome_length = 823786402
inp = open('Sly.1.trimmed.fastq','r')
inl = inp.readline()
while inl:
	if inl.startswith('@'):
		D.append([inl.split(' ')[0][1:]])
	inl = inp.readline()
		
selected_reads = random.sample(D,genome_length * number / 101)

inp = open('Sly.1.trimmed.fastq','r')
inl = inp.readline()
while inl:
	if inl.startswith('@'):
		read = inl.split(' ')[0][1:]
		if read in selected_reads:
			output1.write(inl)
			inl = inp.readline()
			output1.write(inl)
			inl = inp.readline()
			output1.write(inl)
			inl = inp.readline()
			output1.write(inl)
	inl = inp.readline()
	
output1.close()

inp = open('Sly.2.trimmed.fastq','r')
inl = inp.readline()
while inl:
	if inl.startswith('@'):
		read = inl.split(' ')[0][1:]
		if read in selected_reads:
			output2.write(inl)
			inl = inp.readline()
			output2.write(inl)
			inl = inp.readline()
			output2.write(inl)
			inl = inp.readline()
			output2.write(inl)
	inl = inp.readline()
	
output2.close()

			
			
	
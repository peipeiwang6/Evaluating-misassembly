import sys,os

path = sys.argv[1] # path with the html files from tandem_repeats_finder

out = open('Tandem_repeats_whole_genome.txt','w')
out.write('Chr\tstart\tstop\tconsensus_seq\tlength\tcopy_number\n')
for files in os.listdir(path):
	if files.endswith('.txt.html'):
		inp = open(path + files,'r').readlines()
		x = 0
		while x < len(inp):
			inl = inp[x]
			if inl.startswith('Sequence:'):
				chr = inl.strip().split('Sly_')[1]
				x += 1
			if 'Indices:' in inl:
				left = inl.split('Indices: ')[1].split('--')[0]
				right = inl.split('Indices: ')[1].split('--')[1].split()[0]
				x += 1
				inl = inp[x]
				num = inl.split('Copynumber: ')[1].split()[0]
			if 'Consensus pattern ' in inl:
				x += 1
				inl = inp[x]
				seq = inl.strip()
				out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(chr,left,right,seq,len(seq),num))
				out.flush()
			x += 1

out.close()
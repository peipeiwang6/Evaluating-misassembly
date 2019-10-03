from Bio import SeqIO
import os,sys
from collections import Counter
inp=sys.argv[1] # HC, LC, or BC regions
genome=sys.argv[2]
kmer=int(sys.argv[3]) # the K, we tested for 1~6bp

inf=open(inp,"r").readlines()
out_itself=open(inp+str(kmer)+".txt","w")
def get_kmer_list(seq,n):
    result = []
    for i in range(len(seq)-n+1):
        result.append(seq[i:i+n])
    return result

sequence = {}
fasta_sequences = SeqIO.parse(open(genome),'fasta')
for sequences in fasta_sequences:
	seq=sequences.seq
	id= sequences.id
	length=len(sequences.seq)
	sequence[id]=[id,seq,length]

for region in inf:
	reg = region.split("\t")
	typ = reg[0]
	chr = reg[1].split(":")[0]
	start = int(reg[1].split(":")[1].split("-")[0])
	end = int(reg[1].split(":")[1].split("-")[1])
	itself=str(sequence[chr][1][start:end+1]).upper()
	itself_kmer=get_kmer_list(itself,kmer)
	itself_dict=dict(Counter(itself_kmer))
	for m,n in itself_dict.items():
		out_itself.write("%s\t%s\t%s\n"%(region.strip(),m,n))

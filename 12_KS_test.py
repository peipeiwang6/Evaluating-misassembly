'''
input 1: path to your files
input 2: output file name
'''
import sys, os
import pandas as pd
from scipy import stats

path = sys.argv[1]
out = open(sys.argv[2],'w')

for file in os.listdir('./'):
	df = pd.read_csv(path + file, sep='\t', index_col = 0, header = 0)
	df_hc = df[df['RD']=='HC']
	df_lc = df[df['RD']=='LC']
	df_bac = df[df['RD']=='Background']
	x = 1
	while x < len(df.columns):
		try:
			pvalue = stats.kruskal(df_hc.iloc[:,x],df_bac.iloc[:,x],df_lc.iloc[:,x])[1]
		except:
			pvalue = 1
		out.write(df.columns[x] + '\t%s\n'%pvalue)
		out.flush()
		x += 1

out.close()
import sys,os
Chr = sys.argv[1]
old_chr = sys.argv[2]
new_chr = sys.argv[3]
inp = open('%s.fplot'%Chr,'r').readlines()
out = open('%s_corresponding_between_two_assemblies.txt'%Chr,'w')
x = 1
while x < len(inp)/4:
	try:
		tem1 = inp[4*x + 1].strip().split()
		tem2 = inp[4*x + 2].strip().split()
		if float(tem1[2]) >= 95:
			out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(old_chr,tem1[1],tem2[1],new_chr,tem1[0],tem2[0],'forward',int(tem2[1])-int(tem1[1]),tem1[2]))
	except:
		print('Bottom of the file')
	x += 1

inp = open('%s.rplot'%Chr,'r').readlines()
x = 1
while x < len(inp)/4:
	try:
		tem1 = inp[4*x + 1].strip().split()
		tem2 = inp[4*x + 2].strip().split()
		if float(tem1[2]) >= 95:
			out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(old_chr,tem1[1],tem2[1],new_chr,tem1[0],tem2[0],'reverse',int(tem2[1])-int(tem1[1]),tem1[2]))
	except:
		print('Bottom of the file')
	x += 1
	
out.close()
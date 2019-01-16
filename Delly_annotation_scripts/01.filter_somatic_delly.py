#Arg1: delly output file
#Arg2: column number of tumor information in delly output

import sys
print('###filter somatic')
print(sys.argv[1])
in_file=open(sys.argv[1])
tcol=int(sys.argv[2])-1
if tcol == 9:
	ncol =10
elif tcol == 10:
	ncol = 9
else:
	print('Wrong normal column number! exit!')
	sys.exit()
out_file=open(sys.argv[1]+'.somatic','w')
in_line=in_file.readline().strip()
n=m=0
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		tDV=int(in_indi[tcol].split(':')[9])
		tRV=int(in_indi[tcol].split(':')[11])
		nDV=int(in_indi[ncol].split(':')[9])
		nRV=int(in_indi[ncol].split(':')[11])
		if tDV+tRV==0: #filter out
			m=m+1
		elif nDV+nRV>0:  #filter out
			m=m+1
		else:
			out_file.write(in_line+'\n')
			n=n+1
	in_line=in_file.readline().strip()
print(n+m)
print(n)

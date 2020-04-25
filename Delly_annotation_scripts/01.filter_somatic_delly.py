#Arg1: delly output file
#Arg2: column number of tumor information in delly output

#2019-06-11: nDV+nRV>0 -> genotype != '0/0'
#2020-04-26: svtype == 'DEL' and pos2-pos1 < 1000 and nRV > 1 -> filter out
#2020-04-26: (svype != 'DEL' or pos1-pos1 >= 1000) and nDV > 1 -> filter out



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
		chr1=in_indi[0]
		pos1=int(in_indi[1]
		chr2=(in_indi[7].split('CHR2=')[1]).split(';')[0]
		pos2=int((in_indi[7].split('END=')[1]).split(';')[0])
		svtype=(in_indi[7].split('SVTYPE=')[1]).split(';')[0]
		tDV=int(in_indi[tcol].split(':')[9])
		tRV=int(in_indi[tcol].split(':')[11])
		nDV=int(in_indi[ncol].split(':')[9])
		nRV=int(in_indi[ncol].split(':')[11])
		nGT=in_indi[ncol].split(':')[0]
		if tDV+tRV==0: #filter out
			m=m+1
		elif nGT != '0/0': #filter out
			m=m+1
		elif svtype == 'DEL' and pos2-pos1 < 1000 and nRV > 1:
			'blank'
		elif (svype != 'DEL' or pos1-pos1 >= 1000) and nDV > 1:
			'blank'
		else:
			out_file.write(in_line+'\n')
			n=n+1
	in_line=in_file.readline().strip()
print(n+m)
print(n)

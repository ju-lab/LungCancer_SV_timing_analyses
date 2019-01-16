import sys
in_file=open(sys.argv[1]) # absCN
print(sys.argv[1])
out_file=open(sys.argv[1]+'.gen_fi','w')
gender=sys.argv[2]  # XX or XY

in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		abscn=in_indi[4]
		if gender == 'XX' and chr1 == 'Y':
			'blank'
		elif gender == 'XY' and (chr1 == 'X' or chr1 == 'Y') and abscn != 'NA':
			out_file.write('\t'.join(in_indi[0:4])+'\t'+str(float(abscn)/2)+'\n')
		else:
			out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()


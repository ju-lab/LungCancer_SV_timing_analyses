#181231: majCN == 'NA'

import sys,os
print('###'+sys.argv[0])
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.edit','w')
in_line=in_file.readline().strip()
prev_mincn='1'
while in_line:
	in_indi=in_line.split('\t')
	if 'chromosome' in in_indi[0]:
		header_list=['#CHROM','start_pos','end_pos','totCN','majCN','minCN']
		out_file.write('\t'.join(header_list)+'\n')
	else:
		chr1=in_indi[0][1:-1]
		st_pos=in_indi[1]
		ed_pos=in_indi[2]
		totcn=in_indi[9]
		majcn=in_indi[10]
		mincn=in_indi[11]
		if totcn =='1' and majcn =='NA':
			majcn ='1'
			mincn = '0'
		info_list=[chr1, st_pos, ed_pos, totcn, majcn, mincn]
		out_file.write('\t'.join(info_list)+'\n')
	in_line=in_file.readline().strip()

#Arg1: sequenza segments.txt

import sys
def find_near_core(chr1, start_pos, end_pos):
	start_pos=int(start_pos); end_pos=int(end_pos)
	up_dist=100000000
	down_dist=100000000
	for line in core_list:
		core_indi=line.split('\t')
		core_chr=core_indi[0][1:-1]
		core_st=int(core_indi[1])
		core_end=int(core_indi[2])
		core_info='\t'.join([core_indi[9], core_indi[10], core_indi[11]])
		if core_chr==chr1 and core_end <= start_pos and (start_pos - core_end) < up_dist:
			up_dist=start_pos-core_end
			up_info=core_info
		elif core_chr==chr1 and core_st >= end_pos and (core_st - end_pos) < down_dist:
			down_dist=core_st - end_pos
			down_info=core_info
		elif core_chr==chr1 and core_st >= end_pos and (core_st - end_pos) >= down_dist:
			break
		if up_dist < down_dist:
			final_info=up_info
		elif up_dist > down_dist:
			final_info=down_info
	return final_info
			
		
in_file=open(sys.argv[1])  # sequenza segments
print(sys.argv[1])
out_file=open(sys.argv[1]+'.clean','w')
out_file.write('#CHROM\tstart_pos\tend_pos\ttotCN\tmajCN\tminCN\n')

in_line=in_file.readline().strip()  #pass 1st row
in_line=in_file.readline().strip()
core_list=[]
while in_line:
	in_indi=in_line.split('\t')
	chr1=in_indi[0][1:-1]
	st_pos=int(in_indi[1])
	end_pos=int(in_indi[2])
	dist=end_pos-st_pos
	if dist < 1000000:
		'blank'
	else:
		core_list.append(in_line)
	in_line=in_file.readline().strip()

in_file.seek(0)
modified_list=[]
in_line=in_file.readline().strip()
in_line=in_file.readline().strip()
while in_line:
	in_indi=in_line.split('\t')
	chr1=in_indi[0][1:-1]
	st_pos=int(in_indi[1])
	end_pos=int(in_indi[2])
	info_list=[chr1,str(st_pos),str(end_pos),in_indi[9],in_indi[10],in_indi[11]]
	if in_line in core_list:
		modified_list.append('\t'.join(info_list))
	else:
		cn_info = find_near_core(chr1, st_pos, end_pos)
		info_list=[chr1, str(st_pos), str(end_pos), cn_info]
		modified_list.append('\t'.join(info_list))
	in_line=in_file.readline().strip()

prv_chr='';prv_cninfo=''
for line in modified_list:
	line_indi=line.split('\t')
	chr1=line_indi[0]
	st_pos=line_indi[1]
	end_pos=line_indi[2]
	cn_info='\t'.join([line_indi[3], line_indi[4], line_indi[5]])
	if chr1 == prv_chr and cn_info == prv_cninfo:
		final_endpos=end_pos
	elif chr1 != prv_chr or cn_info != prv_cninfo:
		if prv_chr != '':
			info_list=[final_chr, final_stpos, final_endpos, final_cninfo]
			out_file.write('\t'.join(info_list)+'\n')
		final_chr = chr1
		final_stpos=st_pos
		final_endpos=end_pos
		final_cninfo=cn_info
	prv_chr=chr1
	prv_cninfo=cn_info

info_list=[final_chr, final_stpos, final_endpos, final_cninfo]
out_file.write('\t'.join(info_list)+'\n')

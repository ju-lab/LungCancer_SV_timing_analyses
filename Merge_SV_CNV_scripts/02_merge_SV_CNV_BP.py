#arg1:CNV file
#arg2:SV file

import sys
s_ra=20000   # SV CNV bp merge search range
minseg=100


#determine centromere start and end position
cent_file=open('/home/users/sypark/03_Tools/annovar/humandb/hg19_cytoBand.txt')
cent_line=cent_file.readline().strip()
cent_dic={}
cent_group_list=['gvar','stalk','acen']
while cent_line:
	cent_indi=cent_line.split('\t')
	chr1=cent_indi[0][3:]
	pos1=cent_indi[1]
	pos2=cent_indi[2]
	group=cent_indi[4]
	if chr1 not in cent_dic.keys():
		cent_dic[chr1]={}
		cent_dic[chr1]['candi']=[]
		cent_dic[chr1]['start']=''
		cent_dic[chr1]['end']=''
	if group in cent_group_list:
		cent_dic[chr1]['candi'].append(int(pos1))
		cent_dic[chr1]['candi'].append(int(pos2))
	cent_line=cent_file.readline().strip()

for chrom in cent_dic.keys():
	cent_dic[chrom]['start']=min(cent_dic[chrom]['candi'])
	cent_dic[chrom]['end']=max(cent_dic[chrom]['candi'])


#add CNV bp
bp_dic={}
in_file=open(sys.argv[1])
in_line=in_file.readline().strip()
prev_chr='0'
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		if chr1 not in bp_dic.keys():
			bp_dic[chr1]=[]
		st_pos=in_indi[1]
		ed_pos=in_indi[2]
		bp_dic[chr1].append(int(st_pos))
		bp_dic[chr1].append(int(ed_pos))
	in_line=in_file.readline().strip()
#rescue chrY stant and end pos
bp_dic['Y']=[]
bp_dic['Y'].append(13400000)
bp_dic['Y'].append(28800000)

#add SV bp
ref_file=open(sys.argv[2])
ref_line=ref_file.readline().strip()
while ref_line:
	if ref_line[0]=='#':
		'blank'
	else:
		ref_indi=ref_line.split('\t')
		chr1=ref_indi[0]
		pos1=ref_indi[1]
		chr2=ref_indi[2]
		pos2=ref_indi[1]
		for bp in bp_dic[chr1]:
			if abs(int(pos1)-bp) < s_ra:
				del bp_dic[chr1][bp_dic[chr1].index(bp)]
		bp_dic[chr1].append(int(pos1))
		if chr2 != '.':
			for bp in bp_dic[chr2]:
				if abs(int(pos2)-bp) < s_ra:
					del bp_dic[chr2][bp_dic[chr2].index(bp)]
			bp_dic[chr2].append(int(pos2))
	ref_line=ref_file.readline().strip()

chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
	
out_file=open(sys.argv[1]+'.SV_CNV_bp.txt','w')
out_file.write('chrom\tstart.pos\tend.pos\n')
for chrom in chr_list:
	cent_st=cent_dic[chrom]['start']
	cent_ed=cent_dic[chrom]['end']
	prev_bp=sorted(bp_dic[chrom])[0]
	for bp in sorted(bp_dic[chrom])[1:]:
		if bp - prev_bp > minseg:
			if bp < cent_st:
				out_file.write(chrom+'\t'+str(prev_bp+1)+'\t'+str(bp)+'\n')
				prev_bp=bp
			elif prev_bp+1 < cent_st and bp >= cent_st and bp < cent_ed:
				out_file.write(chrom+'\t'+str(prev_bp+1)+'\t'+str(cent_st)+'\n')
				prev_bp=cent_st
			elif prev_bp+1 >= cent_st and prev_bp+1 < cent_ed and bp < cent_ed:
				#out_file.write(chrom+'\t'+str(cent_st+1)+'\t'+str(cent_ed)+'\n')
				prev_bp=cent_ed
			elif prev_bp+1 >= cent_st and prev_bp+1 < cent_ed and bp >= cent_ed:
				out_file.write(chrom+'\t'+str(cent_ed+1)+'\t'+str(bp)+'\n')
				prev_bp=bp
			elif prev_bp+1 < cent_st and bp >= cent_ed:
				out_file.write(chrom+'\t'+str(prev_bp+1)+'\t'+str(cent_st)+'\n')
				out_file.write(chrom+'\t'+str(cent_ed+1)+'\t'+str(bp)+'\n')
				prev_bp=bp
			elif prev_bp+1 >= cent_ed:
				out_file.write(chrom+'\t'+str(prev_bp+1)+'\t'+str(bp)+'\n')
				prev_bp=bp

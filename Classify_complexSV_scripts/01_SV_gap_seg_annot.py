#Arg1: SVfile


import sys,os
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.gap_seg','w')
reffai=sys.argv[2]

overlap_bp=1 # basepair number for allowing overlap
overlap_bp=-1*overlap_bp+1

chr_file=open(reffai)
chr_line=chr_file.readline().strip()
chr_size={}
while chr_line:
	chr_indi=chr_line.split('\t')
	chr_size[chr_indi[0]]=chr_indi[1]
	chr_line=chr_file.readline().strip()

in_line=in_file.readline().strip()
line_list=[]; BP_list=[]
while in_line:
	if in_line[0:4]=='#CHR':
		out_file.write(in_line+'\tgap1\tseg1\tgap2\tseg2\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=in_indi[1]
		chr2=in_indi[2]
		pos2=in_indi[3]
		ter1=in_indi[5].split('to')[0]
		ter2=in_indi[5].split('to')[1]
		cluster_num=in_indi[26]
		line_list.append([in_line, cluster_num])
		BP_list.append([chr1,pos1,ter1])
		if pos2!='.':
			BP_list.append([chr2,pos2,ter2])
	in_line=in_file.readline().strip()


for [in_line, cluster_num] in line_list:
	in_indi=in_line.split('\t')
	chr1=in_indi[0]
	pos1=in_indi[1]
	chr2=in_indi[2]
	pos2=in_indi[3]
	ter1=in_indi[5].split('to')[0]
	ter2=in_indi[5].split('to')[1]
	gap1_list=[]; gap2_list=[]
	seg1_list=[]; seg2_list=[]
	if ter1=='3':
		for [chr3,pos3,ter3] in BP_list:
			if chr1==chr3 and ter3=='5' and  int(pos3)-int(pos1) >= overlap_bp: #check
				gap1_list.append(int(pos3)-int(pos1))
			elif chr1==chr3 and ter3=='5' and int(pos3)-int(pos1)< overlap_bp: #check
				seg1_list.append(int(pos3)-int(pos1))
			else:
				'blank'
		if len(gap1_list)==0:
			gap1=int(chr_size[chr1])-int(pos1)
		else:
			gap1=min(gap1_list)

		if len(seg1_list)==0:
			seg1=int(pos1)
		else:
			seg1=abs(max(seg1_list))
	elif ter1=='5':
		for [chr3, pos3, ter3] in BP_list:
			if chr1==chr3 and ter3=='3' and int(pos1)-int(pos3) >= overlap_bp: #check
				gap1_list.append(int(pos1)-int(pos3))
			elif chr1==chr3 and ter3=='3' and int(pos1)-int(pos3) < overlap_bp: #check
				seg1_list.append(int(pos1)-int(pos3))
			else:
				'blank'

		if len(gap1_list)==0:
			gap1=int(pos1)
		else:
			gap1=min(gap1_list)

		if len(seg1_list)==0:
			seg1=int(chr_size[chr1])-int(pos1)
		else:
			seg1=abs(max(seg1_list))
	if ter2=='3':
		for [chr3,pos3,ter3] in BP_list:
			if chr2==chr3 and ter3=='5' and  int(pos3)-int(pos2)>= overlap_bp:  #check
				gap2_list.append(int(pos3)-int(pos2))
			elif chr2==chr3 and ter3=='5' and int(pos3)-int(pos2)< overlap_bp: #check
				seg2_list.append(int(pos3)-int(pos2))
			else:
				'blank'
		if len(gap2_list)==0:
			gap2=int(chr_size[chr2])-int(pos2)
		else:
			gap2=min(gap2_list)

		if len(seg2_list)==0:
			seg2=int(pos2)
		else:
			seg2=abs(max(seg2_list))
	elif ter2=='5':
		for [chr3, pos3, ter3] in BP_list:
			if chr2==chr3 and ter3=='3' and int(pos2)-int(pos3)>= overlap_bp: #check
				gap2_list.append(int(pos2)-int(pos3))
			elif chr2==chr3 and ter3=='3' and int(pos2)-int(pos3) < overlap_bp: #check
				seg2_list.append(int(pos2)-int(pos3))
			else:
				'blank'

		if len(gap2_list)==0:
			gap2=int(pos2)
		else:
			gap2=min(gap2_list)

		if len(seg2_list)==0:
			seg2=int(chr_size[chr2])-int(pos2)
		else:
			seg2=abs(max(seg2_list))
	if pos2 == '.':
		gap2 ='.'; seg2 = '.'
	out_file.write(in_line+'\t'+str(gap1)+'\t'+str(seg1)+'\t'+str(gap2)+'\t'+str(seg2)+'\n')
		

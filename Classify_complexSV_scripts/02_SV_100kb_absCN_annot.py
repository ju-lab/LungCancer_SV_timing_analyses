import sys
print(sys.argv[1])
in_file=open(sys.argv[1])  #sv file
ref_file=open(sys.argv[2])  #100kbcov.absCN
out_file=open(sys.argv[1]+'.100kbAbsCN','w')
in_line=in_file.readline().strip()
BP_list=[]
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]; pos1=in_indi[1]; chr2=in_indi[2]; pos2=in_indi[3]
		if pos1 != '.':
			BP_list.append(chr1+'\t'+pos1)
		if pos2 != '.':
			BP_list.append(chr2+'\t'+pos2)
	in_line=in_file.readline().strip()

ref_line=ref_file.readline().strip()
absCN_dic={}
while ref_line:
	if ref_line[0]=='#':
		'blank'
	else:
		find_list=[]
		ref_indi=ref_line.split('\t')
		chr1=ref_indi[0]
		pos1=int(ref_indi[1])
		absCN=ref_indi[4]
		for chrompos in BP_list:
			chrom=chrompos.split('\t')[0]; posi=int(chrompos.split('\t')[1])
			if chrom==chr1 and posi >= pos1 and posi < pos1+100000:
				absCN_dic[chrom+'\t'+str(posi)]=absCN
				find_list.append(chrom+'\t'+str(posi))
		BP_list=list(set(BP_list)-set(find_list))
		if len(BP_list) == 0:
			break
	ref_line=ref_file.readline().strip()

in_file.seek(0)
in_line=in_file.readline().strip()
while in_line:
	in_indi=in_line.split('\t')
	info_list=in_indi
	if in_line[0:4]=='#CHR':
		header_list=['100kbAbsCN1','100kbAbsCN2']
		out_file.write('\t'.join(info_list)+'\t'+'\t'.join(header_list)+'\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		chr1=in_indi[0]; pos1=in_indi[1]; chr2=in_indi[2]; pos2=in_indi[3]
		kb1=absCN_dic[chr1+'\t'+pos1]
		if pos2 !='.':
			kb2=absCN_dic[chr2+'\t'+pos2]
		else:
			kb2='.'
		out_file.write('\t'.join(info_list)+'\t'+kb1+'\t'+kb2+'\n')
	in_line=in_file.readline().strip()


import sys,numpy
print(sys.argv[1])
in_file=open(sys.argv[1])  # 100kbcov.absCN
out_file=open(sys.argv[1]+'.chrCN','w')
in_line=in_file.readline().strip()
chrCN_dic={}
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		absCN=in_indi[4]
		if chr1 not in chrCN_dic.keys():
			chrCN_dic[chr1]=[]
		if absCN!='NA':
			chrCN_dic[chr1].append(float(absCN))
	in_line=in_file.readline().strip()

chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

for chrom in chr_list:
	if chrom in chrCN_dic.keys():
		absCN_list=chrCN_dic[chrom]
		out_file.write(chrom+'\t'+str(numpy.median(absCN_list))+'\n')



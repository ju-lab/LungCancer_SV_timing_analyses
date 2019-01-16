#edited 171212
import sys,os
in_file=open(sys.argv[1])
print(sys.argv[1])
dist=5000000
out_file=open(sys.argv[1]+'.cluster','w')

input_lines=[]
in_line=in_file.readline().strip()
n=0
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		input_lines.append(in_line)
		n=n+1
	in_line=in_file.readline().strip()
input_line_no=n


def expand_chain(input_clust):
	global output_clust
	output_clust=[]
	for line1 in input_clust:
		output_clust.append(line1)
		in_indi=line1.split('\t')
		in_chr1=in_indi[0];in_pos1=in_indi[1];in_chr2=in_indi[2];in_pos2=in_indi[3]
		for line2 in input_lines:
			ref_indi=line2.split('\t')
			ref_chr1=ref_indi[0];ref_pos1=ref_indi[1];ref_chr2=ref_indi[2];ref_pos2=ref_indi[3]
			if in_chr1==ref_chr1:
				if abs(int(in_pos1)-int(ref_pos1))<dist:
					output_clust.append(line2)
			if in_chr2==ref_chr1 and in_chr2!='.':
				if abs(int(in_pos2)-int(ref_pos1))<dist:
					output_clust.append(line2)
			if in_chr1==ref_chr2:
				if abs(int(in_pos1)-int(ref_pos2))<dist:
					output_clust.append(line2)
			if in_chr2==ref_chr2 and in_chr2 !='.':
				if abs(int(in_pos2)-int(ref_pos2))<dist:
					output_clust.append(line2)
	output_clust=list(set(output_clust))

		

clust_list=[]
m=0
for a in input_lines:
	clust=[]
	clust.append(a)
	clust_no=len(clust)
	expand_chain(clust)
	while len(output_clust)>clust_no:
		clust_no=len(output_clust)
		expand_chain(output_clust)
	output_clust=sorted(output_clust)
	if output_clust in clust_list:
		'blank'
	else:
		clust_list.append(output_clust)
		m=m+len(output_clust)

output_line_no=m


if input_line_no != output_line_no:
	print('Error: No. of input lines != No. of output lines')

n=0
for c in clust_list:
	n=n+1
	out_file.write('#cluster'+str(n)+'\n')
	for d in c:
		out_file.write(d+'\n')

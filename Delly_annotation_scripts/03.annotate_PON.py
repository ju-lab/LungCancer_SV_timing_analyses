#argv1: /path/to/input
#argv2: /path/to/normal_panel
#argv3: ref.fa.fai

import sys,gzip
print('###'+sys.argv[0])
print(sys.argv[1])
dist=200
bin_size=1000000
size_dic={}
size_file=open(sys.argv[3])
size_line=size_file.readline().strip()
while size_line:
	size_indi=size_line.split('\t')
	size_dic[size_indi[0]]=int(size_indi[1])
	size_line=size_file.readline().strip()

print('making dic...')
ref_dic={}
for chrom in size_dic.keys():
	chrom_size=size_dic[chrom]
	ref_dic[chrom]={}
	for i in range(0, chrom_size/bin_size+1):
		ref_dic[chrom][i]=[]

if sys.argv[2][-3:]=='.gz':
	ref_file=gzip.open(sys.argv[2])
else:
	ref_file=open(sys.argv[2])
ref_line=ref_file.readline().strip()
while ref_line:
	if ref_line[0]=='#':
		'blank'
	else:
		ref_indi=ref_line.split('\t')
		chr1=ref_indi[0]
		pos1=int(ref_indi[1])
		chr2=(ref_indi[7].split(';CHR2=')[1]).split(';')[0]
		pos2=int((ref_indi[7].split(';END=')[1]).split(';')[0])
		ctinfo=(ref_indi[7].split(';CT=')[1]).split(';')[0]
		dv=int(ref_indi[9].split(':')[9])
		rv=int(ref_indi[9].split(':')[11])
		s_name=ref_indi[10]
		ref_dic[chr1][pos1/bin_size].append([chr1,pos1,ctinfo,chr2,pos2,dv+rv,s_name])
	ref_line=ref_file.readline().strip()

print('annotating..')
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.pon','w')
in_line=in_file.readline().strip()
while in_line:
	if in_line[0:4]=='#CHR':
		out_file.write(in_line+'\t'+'PON_SV;sampleN\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		sv_c=0;sample_list=[]
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=int(in_indi[1])
		svtype=in_indi[2][0:3]
		chr2=(in_indi[7].split(';CHR2=')[1]).split(';')[0]
		pos2=int((in_indi[7].split(';END=')[1]).split(';')[0])
		ctinfo=(in_indi[7].split(';CT=')[1]).split(';')[0]
		
		if svtype=='DEL':
			distco=min(dist, round(abs(pos1-pos2)/float(4),0))
		else:
			distco=dist
		if pos1%bin_size < distco and pos1 > bin_size:
			target_list= ref_dic[chr1][(pos1/bin_size)-1]+ref_dic[chr1][pos1/bin_size]
		elif (bin_size - pos1%bin_size) < distco and pos1 < size_dic[chr1]-bin_size:
			target_list= ref_dic[chr1][pos1/bin_size]+ref_dic[chr1][(pos1/bin_size)+1]
		else:
			target_list= ref_dic[chr1][pos1/bin_size]
		for [ref_chr1,ref_pos1,ref_ctinfo,ref_chr2,ref_pos2,rc,s_name] in target_list:
			if chr1==ref_chr1:
				if ctinfo == ref_ctinfo:
					if chr2==ref_chr2:
						if abs(pos1-ref_pos1) < distco and abs(pos2-ref_pos2) < distco:
							sv_c += rc
							sample_list.append(s_name)
		sample_list=list(set(sample_list))
		out_file.write(in_line+'\t'+str(sv_c)+';'+str(len(sample_list))+'\n')
	in_line=in_file.readline().strip()

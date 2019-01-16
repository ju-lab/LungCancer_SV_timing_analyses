#Arg1 earlySNVcount (ampseg.earlysnv.txt)
#Arg2 SV clustered
#Arg3 sampleid
#Arg4 outdir


import sys,collections
distco=5000000
print('###'+sys.argv[0])
print(sys.argv[3])
seg_file=open(sys.argv[1]) 
sv_file=open(sys.argv[2])  
out_file=open(sys.argv[4]+'/'+sys.argv[3]+'.cl_time','w')

# make cluster list
sv_line=sv_file.readline().strip()
cln_list=[]
while sv_line:
	if sv_line[0]=='#': 'blank'
	else:
		sv_indi=sv_line.split('\t')
		cln=int(sv_indi[26])  #clusternum
		cln_list.append(cln)
	sv_line=sv_file.readline().strip()
cln_dic=collections.Counter(cln_list)
cln_list=list(set(cln_list))

# classify for cluter and annot BPnum
header='off'
for cl in cln_list:
	cl_bp_list=[]
	cl_CN_list=[]
	sv_file.seek(0)
	sv_line=sv_file.readline().strip()
	while sv_line:
		if sv_line[0]=='#': 'blank'
		else:
			sv_indi=sv_line.split('\t')
			cln=int(sv_indi[26])
			if cln == cl:
				chr1=sv_indi[0];pos1=sv_indi[1]
				chr2=sv_indi[2];pos2=sv_indi[3]
				cl_bp_list.append([chr1,pos1])
				if pos2 !='.':
					cl_bp_list.append([chr2,pos2])
			else:
				'blank'
		sv_line=sv_file.readline().strip()
	
	cl_svn=cln_dic[cl]
	seg_file.seek(0)
	seg_line=seg_file.readline().strip()
	seg_line_list=[]
	while seg_line:
		if seg_line[0:4]=='#chr' and header=='off':
			out_file.write('#cluster\t'+seg_line[1:]+'\tBPnum\tcluster_SVnum\n')
			header='on'
		elif seg_line[0]=='#':
			'blank'
		else:
			seg_indi=seg_line.split('\t')
			seg_chr=seg_indi[0];seg_pos1=int(seg_indi[1]);seg_pos2=int(seg_indi[2])
			seg_chr_num=int(seg_chr.replace('X','23').replace('Y','24'))
			seg_bp=0
			for [chr1,pos1] in cl_bp_list:
				if chr1==seg_chr and int(pos1) >= seg_pos1 and int(pos1) <= seg_pos2:
					seg_bp +=1
				elif chr1==seg_chr and ( abs(int(pos1)-seg_pos1) < distco or abs(int(pos1)-seg_pos2) < distco):
					seg_bp +=1
			if seg_bp > 0:
				info_list=[str(cl),seg_line, str(seg_bp),str(cl_svn)]
				out_file.write('\t'.join(info_list)+'\n')
		seg_line=seg_file.readline().strip()


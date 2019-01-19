#Arg1: SVfile with 100kbAbsCN annot
#Arg2: median chromosomal copy number


import sys, collections,numpy
svnco= 10 # SV number cutoff for evaluation
ampsdco=1.7
balco=0.40  # proportion of balanced BP for defining chromoplexy or blanced chromthripsis
mchrco=3   # number of chromosomes for defining chromoplexy of multichromosomal chromothripsis

sv_file=open(sys.argv[1])  # SV file *.100kbAbsCN
out_file=open(sys.argv[1]+'.complex_class','w')
chr_file=open(sys.argv[2])  # chrCN file *.chrCN
chrCN_dic={}
chr_line=chr_file.readline().strip()
while chr_line:
	chr_indi=chr_line.split('\t')
	chrCN_dic[chr_indi[0]]=chr_indi[1]
	chr_line=chr_file.readline().strip()

cluster_list=[]
sv_line=sv_file.readline().strip()
while sv_line:
	if sv_line[0:4]=='#CHR':
		'blank'
	else:
		sv_indi=sv_line.split('\t')
		clustern=sv_indi[26]
		cluster_list.append(clustern)
	sv_line=sv_file.readline().strip()

clct_dic = collections.Counter(cluster_list)
cluster_list=list(set(cluster_list))
classi_dic={}

for clst in cluster_list:
	if clct_dic[clst] < svnco:
		'blank'
	else:
		gap_list=[]; chr_list=[];validchr_list=[]; amp_base_list=[]
		blbp=0;totbp=0
		classi=''
		sv_file.seek(0)
		sv_line=sv_file.readline().strip()
		while sv_line:
			if sv_line[0]=='#':
				'blank'
			else:
				sv_indi=sv_line.split('\t')
				chr1=sv_indi[0];pos1=sv_indi[1];chr2=sv_indi[2];pos2=sv_indi[3]
				if pos2=='.':
					dist=-1
					chrcn2='.'
				else:
					dist=int(pos2)-int(pos1)
					chrcn2=chrCN_dic[chr2]
				svtype=sv_indi[6]
				clustern=sv_indi[26]
				gap1=sv_indi[27];seg1=sv_indi[28];gap2=sv_indi[29];seg2=sv_indi[30]
				kbcn1=sv_indi[31];kbcn2=sv_indi[32];chrcn1=chrCN_dic[chr1]
				if clustern == clst:
					totbp+=1
					if pos2 !='.':
						totbp+=1
					if gap1 != '.': 
						if svtype == 'DEL' and int(gap1)==int(dist):
							'blank'
						else:
							gap_list.append(int(gap1))
					if gap2 != '.': 
						if svtype=='DEL' and int(gap2) == int(dist):
							'blank'
						else:
							gap_list.append(int(gap2))
					if chr1 != '.': chr_list.append(chr1)
					if chr2 != '.': chr_list.append(chr2)
					if kbcn1 != 'NA' and chrcn1 !='.':
						if float(kbcn1)-float(chrcn1) < 0:
							amp_base_list.append(0)
						else:
							amp_base_list.append(float(kbcn1)-float(chrcn1))
					if kbcn2 != 'NA' and chrcn2 !='.': 
						if float(kbcn2)-float(chrcn2) < 0:
							amp_base_list.append(0)
						else:
							amp_base_list.append(float(kbcn2)-float(chrcn2))
			sv_line=sv_file.readline().strip()
		gap_dic=collections.Counter(gap_list)
		for gap in gap_dic.keys():
			if gap < 500 and gap_dic[gap]>=2:
				blbp+=(gap_dic[gap]//2)*2
#			chr_dic=collections.Counter(chr_list)
#			for chrom in chr_dic.keys():
#				if chr_dic[chrom] >= valchrco:
#					validchr_list.append(chrom)
		chr_list=list(set(chr_list))
		balpro=blbp/float(totbp)
		ampsd=numpy.std(amp_base_list)
		if balpro >= balco and len(chr_list) >= mchrco:
			classi= 'chromoplexy'
		elif balpro >= balco and len(chr_list) < mchrco:
			classi= 'balanced_chromothripsis'
		elif ampsd >= ampsdco :   
			classi= 'Amplification'
		elif len(chr_list) >= mchrco:
			classi= 'multichromosomal_chromothripsis'
		elif len(chr_list) < mchrco:
			classi= 'classical_chromothripsis'
		classi_dic[clst]=classi

sv_file.seek(0)
sv_line=sv_file.readline().strip()
while sv_line:
	if sv_line[0:4]=='#CHR':
		out_file.write(sv_line+'\tclassification\n')
	elif sv_line[0]=='#':
		out_file.write(sv_line+'\n')
	else:
		sv_indi=sv_line.split('\t')
		clustern=sv_indi[26]
		if clustern not in classi_dic.keys():
			out_file.write(sv_line+'\tNA\n')
		else:
			out_file.write(sv_line+'\t'+classi_dic[clustern]+'\n')
	sv_line=sv_file.readline().strip()
		


#Arg1: sequenza edit
#Arg2: hg19_cytoBand.txt

import sys
print('###'+sys.argv[0])
cutoff = 0.60 #check the cutoff for call

#find centromere position
ref_file=open(sys.argv[2])
ref_line=ref_file.readline().strip()
centromere_dic={}
prv_pqinfo=''
while ref_line:
	ref_indi=ref_line.split('\t')
	chr1=ref_indi[0][3:]
	pqinfo=ref_indi[3][0]
	if pqinfo=='q' and prv_pqinfo=='p':
		centromere_dic[chr1]=ref_indi[1]
	prv_pqinfo=pqinfo
	ref_line=ref_file.readline().strip()

print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.broadAMP','w')
header_list=['#CHROM','start_pos','end_pos','totCN','majCN','minCN','broadAMP']
out_file.write('\t'.join(header_list)+'\n')

def replace_name(string):  
	return (((string.replace('AMPE','AMP')).replace('WCAE','WCA')).replace('AMP','AMPE')).replace('WCA','WCAE')

def extend_seg(line_list):
	out_list=[]
	for n in range(0, len(line_list)):
		if line_list[n][2] != '.':
			out_list.append([line_list[n][0],line_list[n][1],line_list[n][2],line_list[n][3]])
		elif line_list[n][1] < 2:
			out_list.append([line_list[n][0],line_list[n][1],line_list[n][2],line_list[n][3]])
		elif n ==0:
			if line_list[n][0]==line_list[n+1][0] and line_list[n+1][2] !='.':
				out_list.append([line_list[n][0],line_list[n][1],line_list[n+1][2],line_list[n][3]])
			else:
				out_list.append([line_list[n][0],line_list[n][1],line_list[n][2],line_list[n][3]])
		elif n > 0 and n < len(line_list)-1:
			if line_list[n][0]==line_list[n-1][0] and line_list[n-1][2] != '.' and n!=0:
				out_list.append([line_list[n][0],line_list[n][1],line_list[n-1][2],line_list[n][3]])
			elif line_list[n][0]==line_list[n+1][0] and line_list[n+1][2] !='.':
				out_list.append([line_list[n][0],line_list[n][1],line_list[n+1][2],line_list[n][3]])
			else:
				out_list.append([line_list[n][0],line_list[n][1],line_list[n][2],line_list[n][3]])
		elif n == len(line_list)-1:
			if line_list[n][0]==line_list[n-1][0] and line_list[n-1][2] != '.' and n!=0:
				out_list.append([line_list[n][0],line_list[n][1],line_list[n-1][2],line_list[n][3]])
			else:
				out_list.append([line_list[n][0],line_list[n][1],line_list[n][2],line_list[n][3]])
			
	return out_list
			
chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
for chrom in chr_list:
	line_list=[]
	centpos=int(centromere_dic[chrom])
	plow=0;phigh=0;qlow=0;qhigh=0
	in_file.seek(0)
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			chr1=in_indi[0]
			pos1=int(in_indi[1])
			pos2=int(in_indi[2])
			dist=pos2-pos1
			totCN=in_indi[3]
			majCN=in_indi[4]
			if totCN == 'NA' or majCN=='NA':
				in_line=in_file.readline().strip()
				continue
			else:
				totCN=int(totCN)
				majCN=int(majCN)
			if chr1 == chrom:
				if pos2 < centpos:
					if majCN >= 2:
						phigh += dist
					elif majCN < 2:
						plow += dist
				elif pos1 < centpos and pos2 >= centpos:
					if majCN >= 2:
						phigh += centpos-pos1
						qhigh += pos2-centpos
					elif majCN < 2:
						plow += centpos-pos1
						qlow += pos2-centpos
				elif pos1 >= centpos:
					if majCN >= 2:
						qhigh += dist
					elif majCN < 2:
						qlow += dist
		in_line = in_file.readline().strip()

	if plow+phigh == 0:
		pamp = '.'
	elif phigh/float(plow+phigh) >= cutoff:
		pamp = chrom+'pAMP'
	else:
		pamp = '.'
	if qlow+qhigh ==0:
		qamp = '.'
	elif qhigh/float(qlow+qhigh) >= cutoff:
		qamp = chrom+'qAMP'
	else:
		qamp ='.'
	if phigh+plow+qhigh+qlow ==0:
		wca = '.'
	elif pamp != '.' and qamp != '.' and (phigh + qhigh)/float(plow+phigh+qlow+qhigh) >=cutoff:
		wca= chrom+'WCA'
	else:
		wca = '.'
	
	in_file.seek(0)
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			chr1=in_indi[0]
			pos1=int(in_indi[1])
			pos2=int(in_indi[2])
			dist=pos2-pos1
			totCN=in_indi[3]
			majCN=in_indi[4]
			if totCN == 'NA' or majCN=='NA':
				in_line=in_file.readline().strip()
				continue
			else:
				totCN=int(totCN)
				majCN=int(majCN)
			if chr1 == chrom:
				if pos2 < centpos:
					info=pamp
				elif pos1< centpos and pos2 >= centpos:
					if centpos-pos1 > pos2-centpos:
						info=pamp
					elif centpos-pos1 <= pos2-centpos:
						info=qamp
				elif pos1 >= centpos:
					info=qamp
				if wca == chrom+'WCA':
					line_list.append([chr1, majCN, wca, in_line])
			#		out_file.write(in_line+'\t'+wca+'\n')
				else:
					line_list.append([chr1,majCN,info,in_line])
			#		out_file.write(in_line+'\t'+info+'\n')
		in_line=in_file.readline().strip()
	
	while extend_seg(line_list) != line_list:
		line_list=extend_seg(line_list)
	
	for [chr1,majCN,info,line] in line_list:
		out_file.write(line+'\t'+info+'\n')
	

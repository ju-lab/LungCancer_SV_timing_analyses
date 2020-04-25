#Arg1: delly output  # must be sorted. (del 3-5, dup 5-3, inv pos1<pos2, tra chr1<chr2; 1...22,X,Y,MT)
#Arg2: Tumor bam
#Arg3: Normal bam

#2020-04-12 set read limitation as 50000

import sys, pysam
import collections
from scipy.stats import ttest_ind
print('### Find BP')
print(sys.argv[1])
dl_file=open(sys.argv[1])   # delly output
dl_line=dl_file.readline().strip()
tbam_file=pysam.AlignmentFile(sys.argv[2],'rb')  # Cancer bam
nbam_file=pysam.AlignmentFile(sys.argv[3],'rb') #Normal bam
out_file=open(sys.argv[1]+'.BPinfo','w')
fors=700; bacs=100  #cut-off check! bacs must be smaller than fors
r_limit = 50000 # The number of reads that this script searches for analysis

def make_cigartuple(cigarstring):
	cg_num=len(cigarstring)
	lt=''
	cigar_tuple_list=[]
	for n in range(0,cg_num):
		try: lt = lt+str(int(cigarstring[n]))
		except:
			if cigarstring[n]=='M': cigar_tuple_list.append((0,int(lt)))
			elif cigarstring[n]=='I': cigar_tuple_list.append((1,int(lt)))
			elif cigarstring[n]=='D': cigar_tuple_list.append((2,int(lt)))
			elif cigarstring[n]=='N': cigar_tuple_list.append((3,int(lt)))
			elif cigarstring[n]=='S': cigar_tuple_list.append((4,int(lt)))
			elif cigarstring[n]=='H': cigar_tuple_list.append((5,int(lt)))
			elif cigarstring[n]=='P': cigar_tuple_list.append((6,int(lt)))
			elif cigarstring[n]=='=': cigar_tuple_list.append((7,int(lt)))
			elif cigarstring[n]=='X': cigar_tuple_list.append((8,int(lt)))
			elif cigarstring[n]=='B': cigar_tuple_list.append((9,int(lt)))
			lt=''
	return cigar_tuple_list

def find_M_range(cigar):
	m_start=0;m_end=0  # m_start: just before the start, m_end= the exact end
	cigar_list=make_cigartuple(cigar)
	m_count=0
	for (t, n) in cigar_list:
		if t == 0:
			m_count +=1
	
	if m_count ==1:
		for (t,n) in cigar_list:
			if t!=0 and t!=1:
				m_start+=n
			elif t==0:
				m_end=m_start+n
				break
	elif m_count > 1:
		find_m=0;m_length=0
		for (t,n) in cigar_list:
			if find_m==0 and t!=0 and t!=1:
				m_start+=n
			elif find_m >0 and t!=0 and t!=1:
				m_length+=n
			elif t==0:
				find_m+=1
				if find_m < m_count:
					m_length+=n
				elif find_m == m_count:
					m_end=m_start+m_length
					break
	return([m_start, m_end])

def change_chr_to_int(chr1):
	if chr1[0:2]=='GL':
		chr_n=25+float(chr1[2:])
	elif chr1[0:3]=='NC_':
		chr_n=26+(float(chr1[3:])/1000000)
	elif chr1=='hs37d5':
		chr_n = 27
	else:
		chr_n=int(((chr1.replace('X','23')).replace('Y','24')).replace('MT','25'))
	return(chr_n)


def find_SA_reads(chr1,start1, end1, chr2, target_start2, target_end2, bam_file, r_limit):
	saINFO=[]
	reverse_list=['1','3','5','7','9','b','d','f']
	start1=max(start1,1)
	end1=max(end1,1)
	n = 0
	for read in bam_file.fetch(chr1,start1-1,end1):
		n = n+1
		if n > r_limit:
			break
		if read.cigartuples == None or read.is_secondary == True or read.is_supplementary == True or read.is_duplicate == True:
			continue
		if read.has_tag('SA'):
			chr_n=change_chr_to_int(read.reference_name)
			cigar_info=read.cigarstring
			read_size=read.infer_read_length()
			SA_list=str(read.get_tag('SA')).split(';')[:-1]
			for SA_indi in SA_list:
				info_ori=''
				SA_chr=SA_indi.split(',')[0]
				SA_chr_n=change_chr_to_int(SA_chr)
				SA_pos=int(SA_indi.split(',')[1])
				SA_strand=SA_indi.split(',')[2]
				SA_cigar=SA_indi.split(',')[3]
				SA_MQ=SA_indi.split(',')[4]
				if SA_chr==chr2 and int(SA_pos)>=target_start2 and int(SA_pos)<=target_end2: #check
					pri_M_range=find_M_range(cigar_info)
					SA_M_range=find_M_range(SA_cigar)
					pri_len=pri_M_range[1]-pri_M_range[0]
					SA_len=SA_M_range[1]-SA_M_range[0]
					if (hex(int(read.flag))[-2] in reverse_list and SA_strand=='-') or (hex(int(read.flag))[-2] not in reverse_list and SA_strand=='+'): #same_direction
						if pri_M_range[0] <= SA_M_range[0] and pri_M_range[1] >= SA_M_range[1]: continue
						elif SA_M_range[0] <= pri_M_range[0] and SA_M_range[1] >= pri_M_range[1]: continue
						if pri_M_range[1] > SA_M_range[1]:
							MHLEN=SA_M_range[1]-pri_M_range[0]
							bp1=read.reference_start+1
							bp2=SA_pos+SA_len-1
							terminal1="5";terminal2="3"
							if read.reference_name!=SA_chr:
								rearr="TRA"
								if chr_n < SA_chr_n: info_ori='rs'
								elif chr_n > SA_chr_n: info_ori='sr'
							else:
								if bp1<bp2: rearr="DUP"; info_ori='rs'
								elif bp1>bp2: rearr="DEL"; info_ori='sr'

						elif SA_M_range[1] > pri_M_range[1]:
							MHLEN=pri_M_range[1]-SA_M_range[0]
							bp1=read.reference_start+pri_len
							bp2=SA_pos
							terminal1="3"; terminal2="5"
							if read.reference_name!=SA_chr:
								rearr="TRA"
								if chr_n < SA_chr_n: info_ori='rs'
								elif chr_n > SA_chr_n: info_ori='sr'
							else:
								if bp1<bp2: rearr="DEL"; info_ori='rs'
								elif bp1>bp2: rearr="DUP"; info_ori='sr'
						else:
							'blank'
					else:  # opposite direction
						rvs_pri_M_range=[read_size-pri_M_range[1], read_size-pri_M_range[0]]
						if rvs_pri_M_range[0] <= SA_M_range[0] and rvs_pri_M_range[1] >= SA_M_range[1]: continue
						elif SA_M_range[0] <= rvs_pri_M_range[0] and SA_M_range[1] >= rvs_pri_M_range[1]: continue
						if rvs_pri_M_range[1] > SA_M_range[1]:
							MHLEN=SA_M_range[1]-rvs_pri_M_range[0]
							bp1=read.reference_start+pri_len
							bp2=SA_pos+SA_len-1
							terminal1="3";terminal2="3"
							if read.reference_name!=SA_chr:
								rearr="TRA"
								if chr_n < SA_chr_n: info_ori='rs'
								elif chr_n > SA_chr_n: info_ori='sr'
							else:
								rearr="INV"
								if bp1 < bp2: info_ori='rs'
								elif bp1 > bp2:	info_ori='sr'
						elif SA_M_range[1] > rvs_pri_M_range[1]:
							MHLEN=rvs_pri_M_range[1]-SA_M_range[0]
							bp1=read.reference_start+1
							bp2=SA_pos
							terminal1="5";terminal2="5"
							if read.reference_name!=SA_chr:
								rearr="TRA"
								if chr_n < SA_chr_n: info_ori='rs'
								elif chr_n > SA_chr_n: info_ori='sr'
							else:
								rearr="INV"
								if bp1 < bp2: info_ori="rs"
								elif bp1 > bp2:	info_ori="sr"
						else:
							'blank'
				else:
					'blank'	
			if info_ori=='rs':
				rs_info=read.reference_name+':'+str(bp1)+'_'+SA_chr+':'+str(bp2)+'_'+str(MHLEN)+'_'+rearr+'_'+terminal1+'to'+terminal2
				saINFO.append(rs_info)
			elif info_ori=='sr':
				sr_info=SA_chr+':'+str(bp2)+'_'+read.reference_name+':'+str(bp1)+'_'+str(MHLEN)+'_'+rearr+'_'+terminal2+'to'+terminal1
				saINFO.append(sr_info)
		else:
			'blank'		
	return(saINFO)

while dl_line:
	if dl_line[0:2]=='##':
		out_file.write(dl_line+'\n')
	elif dl_line[0:4]=='#CHR':
		out_file.write(dl_line+'\ttBPinfo\tnBPinfo\n')
	elif dl_line[0]=='#':
		out_file.write(dl_line+'\n')
	else:
		dl_indi=dl_line.split('\t')
		chr1=dl_indi[0];pos1=int(dl_indi[1]);chr2=(dl_indi[7].split('CHR2=')[1]).split(';')[0];pos2=int((dl_indi[7].split('END=')[1]).split(';')[0])
		ct1=(dl_indi[7].split('CT=')[1]).split(';')[0][0]
		ct2=(dl_indi[7].split('CT=')[1]).split(';')[0][-1]
		sv_type=dl_indi[2][0:3]
		dist=pos2-pos1
		if sv_type!='TRA' and dist < 0:
			print('Sorting error')
			print(dl_line)
			sys.exit(1)
		if sv_type == 'INV' and ct1=='3' and ct2=='3':
			start1=pos1-bacs
			end1=pos1+fors
			start2=pos2-bacs
			end2=pos2+fors
			if dist >= bacs and dist < fors:
				end1=pos1+dist-bacs
			elif dist < bacs:
				end1=pos1+dist/2
				start2=pos2-dist/2
		elif sv_type == 'INV' and ct1 == '5' and ct2=='5':
			start1=pos1-fors
			end1=pos1+bacs
			start2=pos2-fors
			end2=pos2+bacs
			if dist >= bacs and dist < fors:
				start2=pos2-dist+bacs
			elif dist < bacs:
				end1=pos1+dist/2
				start2=pos2-dist/2
		elif sv_type == 'DEL':
			start1=pos1-bacs
			end1=pos1+fors
			start2=pos2-fors
			end2=pos2+bacs
			if dist < fors:
				end1=pos1+dist/2
				start2=pos2-dist/2
		elif sv_type == 'DUP':
			start1=pos1-fors
			end1=pos1+bacs
			start2=pos2-bacs
			end2=pos2+fors
			if dist < bacs:
				end1=pos1+dist/2
				start2=pos2-dist/2
		elif sv_type == 'INS':
			t_info='NA'; n_info='NA'
			out_file.write(dl_line+'\t'+t_info+'\t'+n_info+'\n')
			dl_line=dl_file.readline().strip()
			continue
		elif sv_type == 'TRA':
			if ct1=='5':
				start1=pos1-fors
				end1=pos1+bacs
			elif ct1=='3':
				start1=pos1-bacs
				end1=pos1+fors
			if ct2=='5':
				start2=pos2-fors
				end2=pos2+bacs
			elif ct2=='3':
				start2=pos2-bacs
				end2=pos2+fors
	
		if sv_type == 'DEL' or sv_type=='DUP':
			target_start1=start1
			target_end1=end1
			target_start2=start2
			target_end2=end2
		elif sv_type=='INV' or sv_type=='TRA':
			if ct1=='5':
				target_start1=pos1-fors
				target_end1=pos1+bacs
			elif ct1=='3':
				target_start1=pos1-bacs
				target_end1=pos1+fors
			if ct2=='5':
				target_start2=pos2-fors
				target_end2=pos2+bacs
			elif ct2=='3':
				target_start2=pos2-bacs
				target_end2=pos2+fors

		t_list1=find_SA_reads(chr1, start1, end1, chr2, target_start2, target_end2,tbam_file, r_limit)
		t_list2=find_SA_reads(chr2, start2, end2, chr1, target_start1, target_end1,tbam_file, r_limit)
		t_list=t_list1+t_list2
		if len(t_list) > 0:
			t_dic=collections.Counter(t_list)
			t_new_list=[]
			for info in t_dic.keys():
				t_new_list.append(info+'('+str(t_dic[info])+')')
			t_info=','.join(t_new_list)
		else:
			t_info='NA'

		n_list1=find_SA_reads(chr1, start1, end1, chr2, target_start2, target_end2, nbam_file, r_limit)
		n_list2=find_SA_reads(chr2, start2, end2, chr1, target_start1, target_end1, nbam_file, r_limit)
		n_list=n_list1+n_list2
		if len(n_list) > 0:
			n_dic=collections.Counter(n_list)
			n_new_list=[]
			for info in n_dic.keys():
				n_new_list.append(info+'('+str(n_dic[info])+')')
			n_info=','.join(n_new_list)
		else:
			n_info='NA'
		out_file.write(dl_line+'\t'+t_info+'\t'+n_info+'\n')
	dl_line=dl_file.readline().strip()

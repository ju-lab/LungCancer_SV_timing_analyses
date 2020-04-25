#Arg1: sv file
#Arg2: tumor bam
#Arg3: normal bam


#2019-06-11 remove condition about supplementary or secondary reads
#2019-10-01 minor error correction
#2020-04-12 add reads number limitation as 50000

import sys,pysam, collections,  itertools
print('### SVvaf new BP')
print(sys.argv[1])
sv_file=open(sys.argv[1]) #CHR1 POS1 CHR2 POS2 MH Terminal SVtype
t_file=pysam.AlignmentFile(sys.argv[2],'rb') #Cancer bam
n_file=pysam.AlignmentFile(sys.argv[3],'rb') #Normal bam
out_file=open(sys.argv[1]+'.SVvaf','w')
fors=700; bacs=5  # search range for pairread counting # You can adjust these values.
iscut=700 # insert size cut for call reference pair # You can adjust this value.
shortDco=500
sc_co=5  # number of bases which is used to differentiate discordant soft-clipping
r_limit = 50000 # The number of reads that this script searches for analysis
sv_line=sv_file.readline().strip()
#Assign the column number starting from 1
c_chr1=15
c_pos1=16
c_chr2=17
c_pos2=18
c_ter=20   # e.g. 3to5, 5to3, 3to3, etc.
c_type=21  # e.g. DEL, TRA, DUP, INV

c_chr1 -=1;c_pos1-=1;c_chr2-=1;c_pos2-=1;c_ter-=1;c_type-=1

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

def estimate_mappedlength(cigarstring):  # this function is same with read.reference_length made for MC tag or SA tag
	cg_num=len(cigarstring)
	lt=''
	current_m=0;current_d=0
	for n in range(0,cg_num):
		try: lt = lt+str(int(cigarstring[n]))
		except:
			if cigarstring[n]=='M': current_m = current_m + int(lt)
			elif cigarstring[n]=='D' and current_m > 0:
				current_d = current_d +int(lt)
			else:'blank'
			lt=''
	return current_m+current_d


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

def find_interCigar_BP(info1, info2, read_size):  # info eg.: 1,162768391,-,43S46M62S
	SA_chr1=info1.split(',')[0]
	SA_chr1n=change_chr_to_int(SA_chr1)
	SA_pos1=int(info1.split(',')[1])
	SA_strand1=info1.split(',')[2]
	SA_cigar1=info1.split(',')[3]
	SA_chr2=info2.split(',')[0]
	SA_chr2n=change_chr_to_int(SA_chr2)
	SA_pos2=int(info2.split(',')[1])
	SA_strand2=info2.split(',')[2]
	SA_cigar2=info2.split(',')[3]
	M_range1=find_M_range(SA_cigar1)
	M_range2=find_M_range(SA_cigar2)
	len1=M_range1[1]-M_range1[0]
	len2=M_range2[1]-M_range2[0]
	if SA_strand1 == SA_strand2: #same_direction
		if M_range1[0] <= M_range2[0] and M_range1[1] >= M_range2[1]: 
			return('overlap')
		elif M_range2[0] <= M_range1[0] and M_range2[1] >= M_range1[1]:
			return('overlap')
		if M_range1[1] > M_range2[1]:
			MHLEN=M_range2[1]-M_range1[0]
			bp1=SA_pos1
			bp2=SA_pos2+len2-1
			terminal1="5";terminal2="3"
			if SA_chr1!=SA_chr2:
				rearr="TRA"
			else:
				if bp1<=bp2: rearr="DUP"
				elif bp1>bp2: rearr="DEL"

		elif M_range2[1] > M_range1[1]:
			MHLEN=M_range1[1]-M_range2[0]
			bp1=SA_pos1+len1-1
			bp2=SA_pos2
			terminal1="3"; terminal2="5"
			if SA_chr1!=SA_chr2:
				rearr="TRA"
			else:
				if bp1<bp2: rearr="DEL"
				elif bp1>=bp2: rearr="DUP"
		else:
			'blank'
	else:  # opposite direction
		rvs_M_range1=[read_size-M_range1[1], read_size-M_range1[0]]
		if rvs_M_range1[0] <= M_range2[0] and rvs_M_range1[1] >= M_range2[1]:
			return('overlap')
		elif M_range2[0] <= rvs_M_range1[0] and M_range2[1] >= rvs_M_range1[1]:
			return('overlap')
		if rvs_M_range1[1] > M_range2[1]:
			MHLEN=M_range2[1]-rvs_M_range1[0]
			bp1=SA_pos1+len1-1
			bp2=SA_pos2+len2-1
			terminal1="3";terminal2="3"
			if SA_chr1!=SA_chr2:
				rearr="TRA"
			else:
				rearr="INV"
		elif M_range2[1] > rvs_M_range1[1]:
			MHLEN=rvs_M_range1[1]-M_range2[0]
			bp1=SA_pos1
			bp2=SA_pos2
			terminal1="5";terminal2="5"
			if SA_chr1!=SA_chr2:
				rearr="TRA"
			else:
				rearr="INV"
		else:
			'blank'
	info=SA_chr1+':'+str(bp1)+';'+SA_chr2+':'+str(bp2)+';'+str(MHLEN)+';'+rearr+';'+terminal1+'to'+terminal2+';'+str(len1)+';'+str(len2)
	rvs_info=SA_chr2+':'+str(bp2)+';'+SA_chr1+':'+str(bp1)+';'+str(MHLEN)+';'+rearr+';'+terminal2+'to'+terminal1+';'+str(len2)+';'+str(len1)
	if SA_chr1n < SA_chr2n: 
		return(info)
	elif SA_chr1n > SA_chr2n:
		return(rvs_info)
	elif SA_chr1n == SA_chr2n:
		if bp1 <= bp2:
			return(info)
		elif bp1 > bp2:
			return(rvs_info)
		

def find_mate_from_SA(read):
	newBP_list=[];neoBP_list=[]
	reverse_list=['1','3','5','7','9','b','d','f']
	cigar_info=read.cigarstring
	read_size=read.infer_read_length()
	SA_list=str(read.get_tag('SA')).split(';')[:-1]
	if hex(int(read.flag))[-2] in reverse_list:
		read_strand='-'
	else:
		read_strand='+'
	read_info=read.reference_name+','+str(read.reference_start+1)+','+read_strand+','+cigar_info
	for SA_indi in SA_list:
		res=find_interCigar_BP(read_info, SA_indi, read_size)
		if res != 'overlap':
			newBP_list.append(res)
	if len(SA_list)>1:
		info_combi=list(itertools.combinations(SA_list,2))
		for (info1,info2) in info_combi:
			res=find_interCigar_BP(info1,info2,read_size)
			if res != 'overlap':
				neoBP_list.append(res)
	return(newBP_list, neoBP_list)

def mate_list_summary(mate_list):
	summary_dic={}
	for mate in mate_list:
		mate_indi=mate.split(';')
		m1=int(mate_indi[5])
		m2=int(mate_indi[6])
		info=';'.join(mate_indi[0:5])
		if info not in summary_dic.keys():
			summary_dic[info]={}
			summary_dic[info]['num']=0
			summary_dic[info]['match1']=[]
			summary_dic[info]['match2']=[]
		summary_dic[info]['num']+=1
		summary_dic[info]['match1'].append(m1)
		summary_dic[info]['match2'].append(m2)
	final_list=[]
	for info in summary_dic.keys():
		m1max=max(summary_dic[info]['match1'])
		m2max=max(summary_dic[info]['match2'])
		freq=summary_dic[info]['num']
		final_list.append(info+';'+str(m1max)+';'+str(m2max)+'('+str(freq)+')')
	return (','.join(final_list))

def find_discordant_reads(chr1, pos1, ter1, chr2, pos2, ter2, pysam_file,sa_seq_list, r_limit):
	pos1=int(pos1); pos2=int(pos2); ter1=int(ter1); ter2=int(ter2)
	new_mate_list=[];neo_mate_list=[]
	if ter1==3: pos1_start=pos1-fors; pos1_end=pos1+bacs
	elif ter1==5: pos1_start=pos1-bacs; pos1_end=pos1+fors
	if ter2==3: pos2_start=pos2-fors; pos2_end=pos2+bacs
	elif ter2==5: pos2_start=pos2-bacs; pos2_end=pos2+fors
	pos1_start=max(pos1_start, 1); pos2_start=max(pos2_start, 1)
	if chr1 == chr2 and ter1==5 and ter2 == 3 and pos1 < pos2:  # exceptional short duplication
		pos1_end=min(pos1_end, pos1+(pos2-pos1)/2)
		pos2_start=max(pos2_start, pos2-(pos2-pos1)/2)
	elif chr1 == chr2 and ter1==3 and ter2 ==5 and pos2 < pos1:
		pos2_end=min(pos2_end, pos2+(pos1-pos2)/2)
		pos1_start=max(pos1_start, pos1-(pos1-pos2)/2)
	pair_true_list=[];sp_true_list=[];sa_true_list=[]
	pair_ref_list=[]; jx_ref_list=[]
	n=0
	for read in pysam_file.fetch(chr1, pos1_start-1, pos1_end):
		n= n+1
		if n > r_limit :
			break
		if read.is_unmapped == True or read.is_paired == False or read.mate_is_unmapped == True or read.is_duplicate == True: continue
		if read.has_tag('SA')== True and ((ter1==3 and (read.cigartuples[-1][0]==4 or read.cigartuples[-1][0]==5)) or (ter1==5 and (read.cigartuples[0][0] ==4 or read.cigartuples[0][0] ==5))):
			SA_list=str(read.get_tag('SA')).split(';')[:-1]
			if read.is_reverse == True: PA_strand='+'
			elif read.is_reverse == False: PA_strand='-'
			SA_BP_candi=[]
			for SA_indi in SA_list:
				SA_chr=SA_indi.split(',')[0]
				SA_pos=int(SA_indi.split(',')[1])
				SA_strand=SA_indi.split(',')[2]
				SA_cigar=SA_indi.split(',')[3]
				SA_cigartuples=make_cigartuple(SA_cigar)
				SA_MQ=SA_indi.split(',')[4]
				current_m=0; current_d=0
				for cigar in SA_cigartuples:
					if cigar[0]==0:
						current_m=current_m+cigar[1]
					elif cigar[0]==2 and current_m > 0:
						current_d=current_d+cigar[1]
					elif (cigar[0]==4 or cigar[0]==5) and current_m > 0:
						break
				if ((SA_cigartuples[0][0]==4 or SA_cigartuples[0][0]==5) and SA_chr == chr2 and abs(SA_pos-pos2) <=1) or ((SA_cigartuples[-1][0]==4 or SA_cigartuples[-1][0]==5) and SA_chr == chr2 and abs(SA_pos + current_m + current_d-1-pos2) <=1):
					sa_true_list.append(read.query_name)
					sp_true_list.append(read.query_name)
					if ter1==3 and read.cigartuples[-1][0]==4:
						sc_seq=read.query_sequence[read.cigartuples[-1][1]*(-1): read.cigartuples[-1][1]*(-1)+sc_co]
						sa_seq_list.append(sc_seq)
						sa_res=find_mate_from_SA(read)	
						new_mate_list=new_mate_list+sa_res[0]
						neo_mate_list=neo_mate_list+sa_res[1]
					elif ter1==5 and read.cigartuples[0][0]==4:
						sc_seq=read.query_sequence[read.cigartuples[0][1]-1-sc_co+1:read.cigartuples[0][1]-1+1]
						sa_seq_list.append(sc_seq)
						sa_res=find_mate_from_SA(read)	
						new_mate_list=new_mate_list+sa_res[0]
						neo_mate_list=neo_mate_list+sa_res[1]
		sa_seq_list=list(set(sa_seq_list))
	n=0
	for read in pysam_file.fetch(chr1, pos1_start-1, pos1_end):
		n = n +1
		if n > r_limit:
			break
		if read.is_unmapped == True or read.is_paired == False or read.mate_is_unmapped == True or read.is_duplicate == True: continue
		pair_ref_mode='off';jx_ref_mode='off'
		if ter1==3:
			if read.is_reverse == False and read.mate_is_reverse == True and read.next_reference_name == chr1 and read.reference_start +1 < pos1 and read.reference_start +1 +read.template_length -1 > pos1 and read.template_length >= 0 and read.template_length < iscut: 
				pair_ref_list.append(read.query_name)
				pair_ref_mode='on'
			if read.reference_start + 1 <= pos1 and read.reference_start + 1 + read.reference_length - 1 > pos1 and read.next_reference_name == chr1:
				jx_ref_list.append(read.query_name)
				jx_ref_mode='on'
			if pair_ref_mode == 'off' and read.is_reverse == False and read.next_reference_name == chr2 and read.next_reference_start +1 >= pos2_start and read.next_reference_start +1 < pos2_end:
				if (ter2==3 and read.mate_is_reverse == False) or (ter2==5 and read.mate_is_reverse == True): 
					pair_true_list.append(read.query_name)
					if read.has_tag('SA') == True and read.query_name not in sa_true_list and read.query_name not in pair_ref_list:
						sa_res=find_mate_from_SA(read)	
						new_mate_list=new_mate_list+sa_res[0]
						neo_mate_list=neo_mate_list+sa_res[1]
						
			if len(sa_seq_list) > 0:
				if pos1 - (read.reference_start +1) +1 == read.reference_length:
					if read.cigartuples[-1][0]==4 and read.cigartuples[-1][1] >= sc_co:
						sc_seq=read.query_sequence[read.cigartuples[-1][1]*(-1): read.cigartuples[-1][1]*(-1)+sc_co]
						if sc_seq in sa_seq_list:
							sp_true_list.append(read.query_name)
		elif ter1==5:
			if read.is_reverse == True and read.mate_is_reverse == False and read.next_reference_name==chr1 and read.reference_start +1 + read.reference_length -1 >= pos1 and read.reference_start + 1 + read.reference_length -1 +read.template_length + 1 < pos1 and read.template_length < 0 and read.template_length*(-1) < iscut:  # in this situation read.template_length is negative value
				pair_ref_list.append(read.query_name)
				pair_ref_mode='on'
			if read.reference_start + 1 < pos1 and read.reference_start + 1 + read.reference_length - 1 >= pos1: 
				jx_ref_list.append(read.query_name)
				jx_ref_mode='on'
			if pair_ref_mode=='off' and read.is_reverse == True and read.next_reference_name == chr2 and read.next_reference_start +1 >= pos2_start and read.next_reference_start +1 < pos2_end:
				if (ter2==3 and read.mate_is_reverse == False) or (ter2==5 and read.mate_is_reverse == True):
					pair_true_list.append(read.query_name)
					if read.has_tag('SA') == True and read.query_name not in sa_true_list and read.query_name not in pair_ref_list:
						sa_res=find_mate_from_SA(read)	
						new_mate_list=new_mate_list+sa_res[0]
						neo_mate_list=neo_mate_list+sa_res[1]

			if len(sa_seq_list) > 0:
				if read.reference_start + 1 == pos1:
					if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= sc_co:
						sc_seq=read.query_sequence[read.cigartuples[0][1]-1-sc_co+1:read.cigartuples[0][1]-1+1]
						if sc_seq in sa_seq_list:
							sp_true_list.append(read.query_name)
	sa_true_list=list(set(sa_true_list))
	pair_ref_list=list(set(pair_ref_list))
	jx_ref_list=list(set(jx_ref_list) & set(pair_ref_list))
	all_ref_list=list(set(pair_ref_list+jx_ref_list)-set(sa_true_list))
	pair_true_list=list(set(pair_true_list)-set(all_ref_list))
	sp_true_list=list(set(sp_true_list))
	all_true_list=list(set(pair_true_list+sp_true_list+sa_true_list))
	if len(new_mate_list)==0:
		new_mate='NA'
	else:
		new_mate=mate_list_summary(new_mate_list)
	if len(neo_mate_list)==0:
		neo_mate='NA'
	else:
		neo_mate=mate_list_summary(neo_mate_list)
	return([pair_true_list, sp_true_list, sa_true_list, pair_ref_list, jx_ref_list, all_ref_list, sa_seq_list, new_mate,neo_mate])


def calc_final_count(chr1, pos1, ter1, chr2, pos2, ter2, t_bam, n_bam, r_limit):
	pos1=int(pos1); pos2=int(pos2); ter1=int(ter1); ter2=int(ter2)
	a1=0; as1=0; asa1=0;r1=0;rj1=0; r2=0; rj2=0; na1=0; nsa1=0
	normal_split1='off';normal_split2='off'
	sa_seq_list=[]
	t1_list=find_discordant_reads(chr1,pos1,ter1,chr2,pos2,ter2,t_bam,sa_seq_list, r_limit)
	sa_seq_list=t1_list[6]
	new_mate1=t1_list[7]
	neo_mate1=t1_list[8]
	n1_list=find_discordant_reads(chr1,pos1,ter1,chr2,pos2,ter2,n_bam,sa_seq_list, r_limit)
	sa_seq_list=[]
	t2_list=find_discordant_reads(chr2,pos2,ter2,chr1,pos1,ter1,t_bam,sa_seq_list, r_limit)
	sa_seq_list=t2_list[6]
	new_mate2=t2_list[7]
	neo_mate2=t2_list[8]
	n2_list=find_discordant_reads(chr2,pos2,ter2,chr1,pos1,ter1,n_bam,sa_seq_list, r_limit)

	n1_pair_list=n1_list[0]
	n1_sp_list=n1_list[1]
	n1_sa_list=n1_list[2]
	n2_pair_list=n2_list[0]
	n2_sp_list=n2_list[1]
	n2_sa_list=n2_list[2]
	t1_pair_list=t1_list[0]
	t1_sp_list=t1_list[1]
	t1_sa_list=t1_list[2]
	t1_rj_list=t1_list[4]
	t1_rt_list=t1_list[5]
	t2_pair_list=t2_list[0]
	t2_sp_list=t2_list[1]
	t2_sa_list=t2_list[2]
	t2_rj_list=t2_list[4]
	t2_rt_list=t2_list[5]

	t1_total_list=list(set(t1_pair_list+t1_sp_list+t1_sa_list))
	t1_sp_list=list(set(t1_sp_list+t1_sa_list))
	t2_total_list=list(set(t2_pair_list+t2_sp_list+t2_sa_list))
	t2_sp_list=list(set(t2_sp_list+t2_sa_list))

	t_tot_n=len(list(set(t1_total_list+t2_total_list)))
	t_split_n=len(list(set(t1_sp_list+t2_sp_list)))
	t_sa_n=len(list(set(t1_sa_list+t2_sa_list)))
	t1_reftot_n=len(t1_rt_list)
	t1_refjx_n=len(t1_rj_list)
	t2_reftot_n=len(t2_rt_list)
	t2_refjx_n=len(t2_rj_list)
	n_tot_n=len(list(set(n1_pair_list+n1_sa_list+n2_pair_list+n2_sa_list)))
	n_sa_n=len(list(set(n1_sa_list+n2_sa_list)))

	return([t_tot_n, t_split_n,t_sa_n, t1_reftot_n, t1_refjx_n, t2_reftot_n, t2_refjx_n, n_tot_n, n_sa_n, new_mate1, neo_mate1, new_mate2, neo_mate2])


def count_frag_num(chr1, pos1, pysam_file, r_limit):
	pos1=int(pos1);total_frag_list=[]
	n=0
	for read in pysam_file.fetch(chr1, pos1-1, pos1):
		n = n+1
		if n > r_limit:
			break
		if read.is_unmapped == True or read.is_paired == False or read.mate_is_unmapped == True or read.is_duplicate == True: continue
		total_frag_list.append(read.query_name)
	total_frag_list=list(set(total_frag_list))
	return len(total_frag_list)



while sv_line:
	if sv_line[0:4]=='#CHR':
		out_file.write(sv_line+'\tTumor_Ref1;Ref2;AllDiscordantFragments;SplitFragments;SATagFragments;Vaf1;Vaf2\tPairNormal_AllFragments;SAFragments;FragCount1;FragCount2\tnew_mate1\tneo_mate1\tnew_mate2\tneo_mate2\n')
	elif sv_line[0]=='#':
		out_file.write(sv_line+'\n')
	else:
		shortDstatus='off'
		sv_indi=sv_line.split('\t')
		chr1=sv_indi[c_chr1]; pos1=sv_indi[c_pos1]; chr2=sv_indi[c_chr2]; pos2=sv_indi[c_pos2]
		svtype=sv_indi[c_type]
		if pos2 == '.' or svtype == 'INS': 
			print_list=[sv_line, '.','.','.','.','.','.']
		else:
			ter1=sv_indi[c_ter].split('to')[0]; ter2=sv_indi[c_ter].split('to')[1]
			res = calc_final_count(chr1, pos1, ter1, chr2, pos2, ter2, t_file, n_file, r_limit)
			adf=res[0];sf=res[1];saf=res[2];ref1=res[3];rj1=res[4];ref2=res[5];rj2=res[6];na1=res[7];nsa1=res[8]
			new_mate1=res[9]; neo_mate1=res[10]; new_mate2=res[11]; neo_mate2=res[12]
			pnfc1=count_frag_num(chr1, pos1, n_file, r_limit)
			pnfc2=count_frag_num(chr2, pos2, n_file, r_limit)
	
			if chr1 == chr2 and ter1 == '3' and ter2 == '5' and abs(int(pos2)-int(pos1)) < shortDco:
				adf=sf
				ref1=rj1
				ref2=rj2
			elif adf == sf:
				ref1=rj1
				ref2=rj2
			if adf+ref1 == 0:
				vaf1 = 'NA'
			else:
				vaf1=str(round((adf)*100/float(adf+ref1),2))+'%'
			if adf+ref2 == 0:
				vaf2 = 'NA'
			else:
				vaf2=str(round((adf)*100/float(adf+ref2),2))+'%'
			# asr1 and asr2 were not counted in 'else' d/t redundancy with r1, r2
			t_info_list=[str(ref1),str(ref2),str(adf),str(sf),str(saf),vaf1,vaf2]
			n_info_list=[str(na1),str(nsa1),str(pnfc1), str(pnfc2)]
			print_list=[sv_line, ';'.join(t_info_list), ';'.join(n_info_list), new_mate1, neo_mate1, new_mate2, neo_mate2]
		out_file.write('\t'.join(print_list)+'\n')
	sv_line=sv_file.readline().strip()

#Arg1: sv file
#Arg2: tumor bam
#Arg3: normal bam


import sys,pysam
print('### find same clipping')
print(sys.argv[1])
sv_file=open(sys.argv[1]) #CHR1 POS1 CHR2 POS2 MH Terminal SVtype
t_file=pysam.AlignmentFile(sys.argv[2],'rb') #Cancer bam
n_file=pysam.AlignmentFile(sys.argv[3],'rb') #Normal bam
out_file=open(sys.argv[1]+'.pnsc','w')

fors=500;bacs=5
n_sr=10 # normal search range 
sc_co=5
sv_line=sv_file.readline().strip()
#Assign the column number starting from 1
c_chr1=15
c_pos1=16
c_chr2=17
c_pos2=18
c_ter=20   # e.g. 3to5, 5to3, 3to3, etc.
c_type=21  # e.g. DEL, TRA, DUP, INV
c_sa=25

c_chr1 -=1;c_pos1-=1;c_chr2-=1;c_pos2-=1;c_ter-=1;c_type-=1;c_sa-=1


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



def find_pnsc(chr1,pos1,ter1, chr2,pos2,ter2, t_file, n_file):
	sa_seq_list=[];sp_true_list=[]
	pos1=int(pos1);pos2=int(pos2)
	if ter1==3: pos1_start=pos1-fors; pos1_end=pos1+bacs
	elif ter1==5: pos1_start=pos1-bacs; pos1_end=pos1+fors
	if ter2==3: pos2_start=pos2-fors; pos2_end=pos2+bacs
	elif ter2==5: pos2_start=pos2-bacs; pos2_end=pos2+fors
	pos1_start=max(pos1_start, 1); pos2_start=max(pos2_start, 1)
	if chr1 == chr2 and ter1==5 and ter2 == 3 and pos1 < pos2:   # exceptional short duplication
		pos1_end=min(pos1_end, pos1+(pos2-pos1)/2)
		pos2_start=max(pos2_start, pos2-(pos2-pos1)/2)
	elif chr1 == chr2 and ter1==3 and ter2 ==5 and pos2 < pos1:
		pos2_end=min(pos2_end, pos2+(pos1-pos2)/2)
		pos1_start=max(pos1_start, pos1-(pos1-pos2)/2)
	for read in t_file.fetch(chr1, pos1_start-1, pos1_end):
		if read.is_unmapped == True or read.is_paired == False or read.mate_is_unmapped == True or read.is_secondary == True or read.is_supplementary == True or read.is_duplicate == True: continue
		if read.has_tag('SA')== True:
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
					if ter1==3 and read.cigartuples[-1][0]==4 and read.cigartuples[-1][1] >= sc_co:
						sc_seq=read.query_sequence[read.cigartuples[-1][1]*(-1): read.cigartuples[-1][1]*(-1)+sc_co]
						sa_seq_list.append(sc_seq)
					elif ter1==5 and read.cigartuples[0][0]==4 and read.cigartuples[0][1] >= sc_co:
						sc_seq=read.query_sequence[read.cigartuples[0][1]-1-sc_co+1:read.cigartuples[0][1]-1+1]
						sa_seq_list.append(sc_seq)
		if ter1==3:
			if read.is_reverse == False and read.next_reference_name == chr2 and read.next_reference_start +1 >= pos2_start and read.next_reference_start +1 < pos2_end:
				if (ter2==3 and read.mate_is_reverse == False) or (ter2==5 and read.mate_is_reverse == True): 
					if read.cigartuples[-1][0]==4 and read.cigartuples[-1][1] >= sc_co:
						sc_seq=read.query_sequence[read.cigartuples[-1][1]*(-1): read.cigartuples[-1][1]*(-1)+sc_co]
						sa_seq_list.append(sc_seq)
		elif ter1==5:
			if read.is_reverse == True and read.next_reference_name == chr2 and read.next_reference_start +1 >= pos2_start and read.next_reference_start +1 < pos2_end:
				if (ter2==3 and read.mate_is_reverse == False) or (ter2==5 and read.mate_is_reverse == True):
					if read.cigartuples[0][0]==4 and read.cigartuples[0][1] >= sc_co:
						sc_seq=read.query_sequence[read.cigartuples[0][1]-1-sc_co+1:read.cigartuples[0][1]-1+1]
						sa_seq_list.append(sc_seq)
	sa_seq_list=list(set(sa_seq_list))
	for read in n_file.fetch(chr1, max(1,pos1-1-n_sr), pos1+n_sr):
		if read.is_unmapped == True or read.is_paired == False or read.mate_is_unmapped == True or read.is_secondary == True or read.is_supplementary == True or read.is_duplicate == True: continue
		if ter1==3:
			if len(sa_seq_list) > 0:
				if read.cigartuples[-1][0]==4 and read.cigartuples[-1][1] >= sc_co:
					sc_seq=read.query_sequence[read.cigartuples[-1][1]*(-1): read.cigartuples[-1][1]*(-1)+sc_co]
					if sc_seq in sa_seq_list:
						sp_true_list.append(read.query_name)
		elif ter1==5:
			if len(sa_seq_list) > 0:
				if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= sc_co:
					sc_seq=read.query_sequence[read.cigartuples[0][1]-1-sc_co+1:read.cigartuples[0][1]-1+1]
					if sc_seq in sa_seq_list:
						sp_true_list.append(read.query_name)
	return len(list(set(sp_true_list)))


while sv_line:
	if sv_line[0:4]=='#CHR':
		out_file.write(sv_line+'\tPairNormalSameClip\n')
	elif sv_line[0]=='#':
		out_file.write(sv_line+'\n')
	else:
		sv_indi=sv_line.split('\t')
		chr1=sv_indi[c_chr1]; pos1=int(sv_indi[c_pos1]); chr2=sv_indi[c_chr2]; pos2=int(sv_indi[c_pos2])
		svtype=sv_indi[c_type]; tsa=sv_indi[c_sa]
		ter1=sv_indi[c_ter].split('to')[0]; ter2=sv_indi[c_ter].split('to')[1]
		if svtype == 'INS':
			pnsc=0
		else:
			ter1=int(ter1);ter2=int(ter2)
			pnsc1=find_pnsc(chr1,pos1,ter1,chr2,pos2,ter2,t_file,n_file)
			pnsc2=find_pnsc(chr2,pos2,ter2,chr1,pos1,ter1,t_file,n_file)
			pnsc=pnsc1+pnsc2
		out_file.write(sv_line+'\t'+str(pnsc)+'\n')
	sv_line=sv_file.readline().strip()

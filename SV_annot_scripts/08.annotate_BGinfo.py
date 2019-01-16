#Arg1 sv_file
#Arg2 tumor_bam
#Arg2 normal_bam


import sys, pysam
print('### Annotate BGinfo')
print(sys.argv[1])
sv_file=open(sys.argv[1])
tbam_file=pysam.AlignmentFile(sys.argv[2],'rb')
nbam_file=pysam.AlignmentFile(sys.argv[3],'rb')
out_file=open(sys.argv[1]+'.bginfo','w')

ser=300 # search range from breakpoint
iscut=1000 # insert size cutoff for normal pair
mate_bin=100000 #mate location binning range
#Assign the column number starting from 1
c_chr1=15
c_pos1=16
c_chr2=17
c_pos2=18
c_ter=20   # e.g. 3to5, 5to3, 3to3, etc.
c_type=21  # e.g. DEL, TRA, DUP, INV

c_chr1-=1;c_pos1-=1;c_chr2-=1;c_pos2-=1;c_ter-=1;c_type-=1

def amount_discordant(chr1, pos1, pysam_file):
	disc_mate_dic={};normal_frag_list=[];total_frag_list=[];disc_frag_list=[];tra_frag_list=[]
	pos1=int(pos1)
	pos1_start=max(pos1-ser, 1)
	pos1_end=pos1+ser
	for read in pysam_file.fetch(chr1, pos1_start-1, pos1_end):
		if read.is_unmapped == True or read.is_paired == False or read.mate_is_unmapped == True or read.is_secondary == True or read.is_supplementary == True or read.is_duplicate == True: continue
		if read.mapping_quality < 1: continue
		total_frag_list.append(read.query_name)
		if read.next_reference_name == chr1 and read.is_reverse == False and read.template_length > 0 and read.template_length < iscut:
			normal_frag_list.append(read.query_name)
			continue
		elif read.next_reference_name == chr1 and read.is_reverse == True and read.template_length < 0 and read.template_length*(-1) < iscut:
			normal_frag_list.append(read.query_name)
			continue
		disc_frag_list.append(read.query_name)
		if read.next_reference_name != chr1:
			tra_frag_list.append(read.query_name)
		if read.next_reference_name not in disc_mate_dic.keys():
			disc_mate_dic[read.next_reference_name]={}
		binn=(read.next_reference_start+1)/mate_bin
		if binn not in disc_mate_dic[read.next_reference_name].keys():
			disc_mate_dic[read.next_reference_name][binn]=[]
		disc_mate_dic[read.next_reference_name][binn].append(read.query_name)
	total_fn=len(list(set(total_frag_list)))
	normalp_fn=len(list(set(normal_frag_list)))
	disc_fn=len(list(set(disc_frag_list)))
	tra_fn=len(list(set(tra_frag_list)))
	disc_chr_n=len(disc_mate_dic)
	disc_bin_n=0;disc_bin_n2=0
	disc_chr_n2_list=[]
	for chrom in disc_mate_dic.keys():
		disc_bin_n += len(disc_mate_dic[chrom])
		if len(disc_mate_dic[chrom]) >= 2 :
			disc_chr_n2_list.append(chrom)
		for eachbin in disc_mate_dic[chrom].keys():
			if len(disc_mate_dic[chrom][eachbin]) >=2:
				disc_bin_n2 +=1
				disc_chr_n2_list.append(chrom)
	disc_chr_n2=len(list(set(disc_chr_n2_list)))
	info_list=[str(total_fn), str(normalp_fn), str(disc_fn), str(disc_chr_n), str(disc_bin_n), str(disc_chr_n2), str(disc_bin_n2)]
	return ';'.join(info_list)


sv_line=sv_file.readline().strip()
while sv_line:
	if sv_line[0:4]=='#CHR':
		out_file.write(sv_line+'\tTumor_BP1_Total;Normal;Discor;Chr;clust;Chr2;clust2\tTumor_BP2_Total;Normal;Discor;Chr;clust;Chr2;clust2\tNormal_BP1_Total;Normal;Discor;Chr;clust;Chr2;clust2\tNormal_BP2_Total;Normal;Discor;Chr;clust;Chr2;clust2\n')
	elif sv_line[0]=='#':
		out_file.write(sv_line+'\n')
	else:
		sv_indi=sv_line.split('\t')
		chr1=sv_indi[c_chr1]
		pos1=sv_indi[c_pos1]
		chr2=sv_indi[c_chr2]
		pos2=sv_indi[c_pos2]
		res1=amount_discordant(chr1,pos1,tbam_file)
		res2=amount_discordant(chr2,pos2,tbam_file)
		res3=amount_discordant(chr1,pos1,nbam_file)
		res4=amount_discordant(chr2,pos2,nbam_file)
		out_list=[sv_line,res1,res2,res3,res4]
		out_file.write('\t'.join(out_list)+'\n')
	sv_line=sv_file.readline().strip()

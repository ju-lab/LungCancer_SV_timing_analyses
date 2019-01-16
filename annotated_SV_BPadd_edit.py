#Arg1 filtered SV
#Arg2 normal_bam

#181218: add exception of finding newSAmate; add filtering of SAits with out close opposite for bridging; add breakpoint when mate depth is under the cutoff (current 3000)
#181219: removal of newSAmate when initial SV is exsit ; oppo_candi_list[0][0] < allfrag/20; oppo_candi, med_mapq < 15; its_depth1 and depth2 < 3000
#181228: when find additional SV (new and neo), length of soft clipping should be more than 30
#181229: when find additional SV (new and neo), the SV should have SA tage number 2 or more.

import sys,pysam
from operator import itemgetter
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.BPedit','w')
n_file=pysam.AlignmentFile(sys.argv[2], 'rb')
header_list=['#CHR1','POS1','CHR2','POS2','mh','terinfo','svtype','precise','class1','conc1','allfrag','spfrag','safrag','ref1','ref2','mapq1','mapq2','tbginfo1','tbginfo2','nbginfo2','nbginfo2','pos_dist1','pos_dist2','dellyid','unknown_bp1','unknown_bp2']
out_file.write('\t'.join(header_list)+'\n')

s_ra=1000 #search range
SA_bp_v=10  # SA breakpoint variation cutoff

initial_list=[]
initial_id_list=[]
output_list=[]
output_id_list=[]


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

## add precise SV to initial_list
n=0
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		dellyid=in_indi[2]
                chr1=in_indi[14]
		pos1=in_indi[15]
		chr2=in_indi[16]
		pos2=in_indi[17]
		mh=in_indi[18]
		terinfo=in_indi[19]
		svtype=in_indi[20]
		if svtype=='INS':
			in_line=in_file.readline().strip()
			continue
		tfinfo=in_indi[26]
		nfinfo=in_indi[27]
		new1=in_indi[28]
		neo1=in_indi[29]
		new2=in_indi[30]
		neo2=in_indi[31]
		tbginfo1=in_indi[32]
		tbginfo2=in_indi[33]
		nbginfo1=in_indi[34]
		nbginfo2=in_indi[35]
		pnsc=in_indi[36]
		mapq1=in_indi[37]
		mapq2=in_indi[38]
		pos_dist1=in_indi[39]
		pos_dist2=in_indi[40]

		#processing
		idx='\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype])
		ter1=terinfo.split('to')[0]
		ter2=terinfo.split('to')[1]
		ref1=(tfinfo.split(';')[0])
		ref2=(tfinfo.split(';')[1])
		allfrag=int(tfinfo.split(';')[2])
		spfrag=int(tfinfo.split(';')[3])
		safrag=int(tfinfo.split(';')[4])
		nbgbint1=int(nbginfo1.split(';')[4])
		nbgbint2=int(nbginfo2.split(';')[4])
		if mapq1 == 'NA':
			med_mapq1=60
			max_mapq1=60
		else:
			med_mapq1=float(mapq1.split(';')[1])
			max_mapq1=float(mapq1.split(';')[2])
		if mapq2 == 'NA':
			med_mapq2=60
			max_mapq2=60
		else:
			med_mapq2=float(mapq2.split(';')[1])
			max_mapq2=float(mapq2.split(';')[2])
		
		# conditions
		check_list=[]
		unknown_bp1='F'
		unknown_bp2='F'
		if nbgbint1 >=10 or max_mapq1 < 10:
			unknown_bp1='T'
		if nbgbint2 >=10 or max_mapq2 < 10:
			unknown_bp2='T'
	
		if mh != '.':
			n+=1
			precise='Y'
			class1='initial'
			conc1='unchanged'
			final_info_list=[chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
			initial_list.append(final_info_list)
			initial_id_list.append(idx)
	in_line=in_file.readline().strip()


#add non-precise SV to initial_list, if there is no duplicate.
in_file.seek(0)
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		dellyid=in_indi[2]
                chr1=in_indi[14]
		pos1=in_indi[15]
		chr2=in_indi[16]
		pos2=in_indi[17]
		mh=in_indi[18]
		terinfo=in_indi[19]
		svtype=in_indi[20]
		if svtype=='INS':
			in_line=in_file.readline().strip()
			continue
		tfinfo=in_indi[26]
		nfinfo=in_indi[27]
		new1=in_indi[28]
		neo1=in_indi[29]
		new2=in_indi[30]
		neo2=in_indi[31]
		tbginfo1=in_indi[32]
		tbginfo2=in_indi[33]
		nbginfo1=in_indi[34]
		nbginfo2=in_indi[35]
		pnsc=in_indi[36]
		mapq1=in_indi[37]
		mapq2=in_indi[38]
		pos_dist1=in_indi[39]
		pos_dist2=in_indi[40]

		#processing
		idx='\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype])
		ter1=terinfo.split('to')[0]
		ter2=terinfo.split('to')[1]
		ref1=(tfinfo.split(';')[0])
		ref2=(tfinfo.split(';')[1])
		allfrag=int(tfinfo.split(';')[2])
		spfrag=int(tfinfo.split(';')[3])
		safrag=int(tfinfo.split(';')[4])
		nbgbint1=int(nbginfo1.split(';')[4])
		nbgbint2=int(nbginfo2.split(';')[4])

		if mh == '.':
			n+=1
			precise='N'
			class1='initial'
			conc1='unchanged'
			final_info_list=[chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
			ini_search = 'NE'
			for ini_info in initial_id_list:
				ini_chr1=ini_info.split('\t')[0]
				ini_pos1=int(ini_info.split('\t')[1])
				ini_chr2=ini_info.split('\t')[2]
				ini_pos2=int(ini_info.split('\t')[3])
				ini_mh=ini_info.split('\t')[4]
				ini_tinfo=ini_info.split('\t')[5]
				ini_stype=ini_info.split('\t')[6]
				if chr1 == ini_chr1 and abs(int(pos1)-ini_pos1) <= SA_bp_v and chr2==ini_chr2 and abs(int(pos2)-ini_pos2) <= SA_bp_v and terinfo==ini_tinfo:
					ini_search='E'
					break
			if ini_search == 'NE':
				initial_list.append(final_info_list)
				initial_id_list.append(idx)
	in_line=in_file.readline().strip()
print(n)


# add new, neo SV to initial list unless it is in unknown region.
in_file.seek(0)
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		dellyid=in_indi[2]
                chr1=in_indi[14]
		pos1=in_indi[15]
		chr2=in_indi[16]
		pos2=in_indi[17]
		mh=in_indi[18]
		terinfo=in_indi[19]
		svtype=in_indi[20]
		if svtype=='INS':
			in_line=in_file.readline().strip()
			continue
		tfinfo=in_indi[26]
		nfinfo=in_indi[27]
		new1=in_indi[28]
		neo1=in_indi[29]
		new2=in_indi[30]
		neo2=in_indi[31]
		tbginfo1=in_indi[32]
		tbginfo2=in_indi[33]
		nbginfo1=in_indi[34]
		nbginfo2=in_indi[35]
		pnsc=in_indi[36]
		mapq1=in_indi[37]
		mapq2=in_indi[38]
		pos_dist1=in_indi[39]
		pos_dist2=in_indi[40]

		#processing
		idx='\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype])
		ter1=terinfo.split('to')[0]
		ter2=terinfo.split('to')[1]
		ref1=(tfinfo.split(';')[0])
		ref2=(tfinfo.split(';')[1])
		allfrag=int(tfinfo.split(';')[2])
		spfrag=int(tfinfo.split(';')[3])
		safrag=int(tfinfo.split(';')[4])
		nbgbint1=int(nbginfo1.split(';')[4])
		nbgbint2=int(nbginfo2.split(';')[4])
		if mapq1 == 'NA':
			med_mapq1=60
			max_mapq1=60
		else:
			med_mapq1=float(mapq1.split(';')[1])
			max_mapq1=float(mapq1.split(';')[2])
		if mapq2 == 'NA':
			med_mapq2=60
			max_mapq2=60
		else:
			med_mapq2=float(mapq2.split(';')[1])
			max_mapq2=float(mapq2.split(';')[2])
		
		# conditions
		check_list=[]
		unknown_bp1='F'
		unknown_bp2='F'
		if nbgbint1 >=10 or max_mapq1 < 10:
			unknown_bp1='T'
		if nbgbint2 >=10 or max_mapq2 < 10:
			unknown_bp2='T'
		if unknown_bp1 == 'F':
			check_list.append([new1,neo1,chr1,pos1,ter1,med_mapq1,chr2,pos2,ter2])
		if unknown_bp2 == 'F':
			check_list.append([new2,neo2,chr2,pos2,ter2,med_mapq2,chr1,pos1,ter1])

		for [new1,neo1,chr1,pos1,ter1,med_mapq1,chr2,pos2,ter2] in check_list:
			oppo_candi_list=[]
			if neo1 == 'NA':
				neo1n=0
			else:
				neo1n=len(neo1.split(','))
			new1sa_max=0;sa_mate_n1=0;tsan1=0   # tsan = target SA number
			if new1 != 'NA':
				tom_candi_list=[] # target other mate list
				osan_list=[]  # other SA number list
				otsan=0 
				for info in new1.split(','):
					target='off'
					san=int((info.split('(')[1]).split(')')[0])
					c1=info.split(':')[0]
					p1=(info.split(';')[0]).split(':')[1]
					c2=(info.split(';')[1]).split(':')[0]
					p2=(info.split(';')[1]).split(':')[1]
					m1=info.split(';')[2]
					s1=info.split(';')[3]
					tinfo1=info.split(';')[4]
					t1=(info.split(';')[4]).split('to')[0]
					t2=(info.split(';')[4]).split('to')[1]
					len_clip1=int(info.split(';')[5])
					len_clip2=int((info.split(';')[6]).split('(')[0])
					p1=int(p1); p2=int(p2)
					pos1=int(pos1);pos2=int(pos2)
					if c1 == chr1 and c2 !=chr1:
						proxy="1"
					elif c1 != chr1 and c2 == chr1:
						proxy="2"
					elif c1 == chr1 and c2 == chr1:
						if abs(pos1-p1) < abs(pos1-p2):
							proxy="1"
						elif abs(pos1-p1) > abs(pos1-p2):
							proxy="2"
					if proxy=='1':
						len_info=len_clip2
					elif proxy=='2':
						len_info=len_clip1
					if c1==chr1 and t1==ter1 and abs(pos1-p1) <= SA_bp_v:
						target='on'
						if c2 == chr2 and t2==ter2 and abs(pos2-p2) <= SA_bp_v:
							tsan1+=san
							final_info_list=[c1,str(p1),c2,str(p2),m1,tinfo1,s1,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
						elif c2 == chr2 and ((t2 == '3' and ter2 == '3' and p2-pos2 < 0 and p2-pos2 >= -1*s_ra) or (t2 == '5' and ter2 == '5' and pos2-p2 < 0 and pos2-p2 >= -1*s_ra)):
							'blank'
						else:
							precise='Y'
							class1='newSAmate('+str(len_clip1)+','+str(len_clip2)+')'
							conc1='unchanged'
							final_info_list=[c1,str(p1),c2,str(p2),m1,tinfo1,s1,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
							tom_candi_list.append([san,len_info,proxy,final_info_list])
					elif c2==chr1 and t2==ter1 and abs(pos1-p2) <= SA_bp_v:
						target='on'
						if c1 == chr2 and t1 == ter2 and abs(p1-pos2) <= SA_bp_v:
							tsan1 += san
							final_info_list=[c1,str(p1),c2,str(p2),m1,tinfo1,s1,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
						elif c1 == chr1 and ((t1 == '3' and ter1 == '3' and p1-pos1 < 0 and p1-pos1 >= -1*s_ra) or (t1 == '5' and ter1 == '5' and pos1-p1 < 0 and pos1-p1 >= -1*s_ra)):
							'blank'
						else:
							precise='Y'
							class1='newSAmate('+str(len_clip1)+','+str(len_clip2)+')'
							conc1='unchanged'
							final_info_list=[c1,str(p1),c2,str(p2),m1,tinfo1,s1,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
							tom_candi_list.append([san,len_info,proxy,final_info_list])
					if target == 'off' and ((proxy=='1' and t1 != ter1) or (proxy == '2' and t2 != ter1)) and san > 0:
						precise='Y'
						class1='oppoSA('+str(len_clip1)+','+str(len_clip2)+')'
						conc1='unchanged'
						final_info_list=[c1,str(p1),c2,str(p2),m1,tinfo1,s1,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
						oppo_candi_list.append([san,len_info,proxy,final_info_list])
				for oppo_candi_list in [oppo_candi_list, tom_candi_list]:

					if len(oppo_candi_list) >0:
						oppo_candi_list.sort(key=itemgetter(0,1), reverse=True)
						selected_oppo=oppo_candi_list[0][3]
						oppo_proxy=oppo_candi_list[0][2]
						if oppo_candi_list[0][0] < allfrag/20 or oppo_candi_list[0][0] < 2 or med_mapq1 < 15 or oppo_candi_list[0][1] <30:
							continue
						if len(oppo_candi_list) >=3 and oppo_candi_list[0][0] ==1:
							if oppo_candi_list[0][2] == '1':
								selected_oppo=selected_oppo[0:25]+['T']
							elif oppo_candi_list[0][2] == '2':
								selected_oppo=selected_oppo[0:24]+['T']+[selected_oppo[25]]
						oppo_chr1=selected_oppo[0]
						oppo_pos1=int(selected_oppo[1])
						oppo_chr2=selected_oppo[2]
						oppo_pos2=int(selected_oppo[3])
						oppo_mh=selected_oppo[4]
						oppo_tinfo=selected_oppo[5]
						oppo_stype=selected_oppo[6]
						if oppo_proxy == '1':
							oppo_depth_array=n_file.count_coverage(oppo_chr2, int(oppo_pos2)-1, int(oppo_pos2), quality_threshold=0)
						elif oppo_proxy == '2':
							oppo_depth_array=n_file.count_coverage(oppo_chr1, int(oppo_pos1)-1, int(oppo_pos1), quality_threshold=0)
						oppo_depth=oppo_depth_array[0][0]+oppo_depth_array[1][0]+oppo_depth_array[2][0]+oppo_depth_array[3][0]
						ini_search = 'NE'
						for ini_info in initial_id_list:
							ini_chr1=ini_info.split('\t')[0]
							ini_pos1=int(ini_info.split('\t')[1])
							ini_chr2=ini_info.split('\t')[2]
							ini_pos2=int(ini_info.split('\t')[3])
							ini_mh=ini_info.split('\t')[4]
							ini_tinfo=ini_info.split('\t')[5]
							ini_stype=ini_info.split('\t')[6]
							if oppo_chr1 == ini_chr1 and abs(oppo_pos1-ini_pos1) <= SA_bp_v and oppo_chr2==ini_chr2 and abs(oppo_pos2-ini_pos2) <= SA_bp_v and oppo_tinfo==ini_tinfo:
								ini_search='E'
								break
						if ini_search == 'NE' and oppo_depth < 3000:   # depth filter check
							initial_list.append(selected_oppo)
							initial_id_list.append('\t'.join(selected_oppo[0:7]))
				if neo1 != 'NA':
					for info in neo1.split(','):
						san=int((info.split('(')[1]).split(')')[0])
						c1=info.split(':')[0]
						p1=int((info.split(';')[0]).split(':')[1])
						c2=(info.split(';')[1]).split(':')[0]
						p2=int((info.split(';')[1]).split(':')[1])
						m1=info.split(';')[2]
						s1=info.split(';')[3]
						tinfo1=info.split(';')[4]
						len_clip1=int(info.split(';')[5])
						len_clip2=int((info.split(';')[6]).split('(')[0])
						precise='Y'
						class1='SAits('+str(len_clip1)+','+str(len_clip2)+')'
						conc1='unchanged'
						its_depth_array1=n_file.count_coverage(c1, int(p1)-1, int(p1), quality_threshold=0)
						its_depth1=its_depth_array1[0][0]+its_depth_array1[1][0]+its_depth_array1[2][0]+its_depth_array1[3][0]
						its_depth_array2=n_file.count_coverage(c2, int(p2)-1, int(p2), quality_threshold=0)
						its_depth2=its_depth_array2[0][0]+its_depth_array2[1][0]+its_depth_array2[2][0]+its_depth_array2[3][0]
						ini_search = 'NE'
						for ini_info in initial_id_list:
							ini_chr1=ini_info.split('\t')[0]
							ini_pos1=int(ini_info.split('\t')[1])
							ini_chr2=ini_info.split('\t')[2]
							ini_pos2=int(ini_info.split('\t')[3])
							ini_mh=ini_info.split('\t')[4]
							ini_tinfo=ini_info.split('\t')[5]
							ini_stype=ini_info.split('\t')[6]
							if c1 == ini_chr1 and abs(p1-ini_pos1) <= SA_bp_v and c2==ini_chr2 and abs(p2-ini_pos2) <= SA_bp_v and tinfo1==ini_tinfo:
								ini_search='E'
								break
						if ini_search == 'NE' and its_depth1 < 3000 and its_depth2 < 3000 and min(len_clip1, len_clip2) >=30 and san>=2:
							final_info_list=[c1,str(p1),c2,str(p2),m1,tinfo1,s1,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
							initial_list.append(final_info_list)
							initial_id_list.append('\t'.join([c1,str(p1),c2,str(p2),m1,tinfo1,s1]))
	in_line=in_file.readline().strip()



# marking removal to non-precise SV, if there are replacable precise SV in both BPs. Supporting read count of removed SV is merged in to precise SV which replace it. Modified SVs are added to output_list
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in initial_list:
	candidate_info1=''
	candidate_info2=''
	c_dist1=10000
	c_dist2=10000
	i_dist1=10000
	i_dist2=10000
	nbgbint1=int(nbginfo1.split(';')[-1])
	nbgbint2=int(nbginfo2.split(';')[-1])
	if mapq1 == 'NA':
		med_mapq1=60
		max_mapq1=60
	else:
		med_mapq1=float(mapq1.split(';')[1])
		max_mapq1=float(mapq1.split(';')[2])
	if mapq2 == 'NA':
		med_mapq2=60
		max_mapq2=60
	else:
		med_mapq2=float(mapq2.split(';')[1])
		max_mapq2=float(mapq2.split(';')[2])
	# definition of unknown BP
	unknown_bp1='F'
	unknown_bp2='F'
	if nbgbint1 >=10 or max_mapq1 < 10:
		unknown_bp1='T'
	if nbgbint2 >=10 or max_mapq2 < 10:
		unknown_bp2='T'

	if precise=='N':
		for [chr1_2,pos1_2,chr2_2,pos2_2,mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,allfrag_2,spfrag_2,safrag_2,ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2] in initial_list:
			if precise_2=='Y':
				if chr1==chr1_2 and terinfo[0]=='3' and terinfo_2[0]=='3':
					i_dist1=int(pos1_2)-int(pos1) 
				elif chr1==chr1_2 and terinfo[0]=='5' and terinfo_2[0]=='5':
					i_dist1=int(pos1)-int(pos1_2)
				if i_dist1 < s_ra and i_dist1 > -1*SA_bp_v and i_dist1  < c_dist1:
					candidate_info1=[chr1_2,int(pos1_2),chr2_2,int(pos2_2),mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,str(allfrag_2),str(spfrag_2),str(safrag_2),ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2]
					c_dist1=i_dist1
				if chr1==chr2_2 and terinfo[0]=='3' and terinfo_2[-1]=='3':
					i_dist1=int(pos2_2)-int(pos1)
				elif chr1==chr2_2 and terinfo[0]=='5' and terinfo_2[-1]=='5':
					i_dist1=int(pos1)-int(pos2_2)
				if i_dist1 < s_ra and i_dist1 > -1*SA_bp_v and i_dist1  < c_dist1:
					candidate_info1=[chr1_2,int(pos1_2),chr2_2,int(pos2_2),mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,str(allfrag_2),str(spfrag_2),str(safrag_2),ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2]
					c_dist1=i_dist1
	
				if chr2==chr1_2 and terinfo[-1]=='3' and terinfo_2[0]=='3':
					i_dist2=int(pos1_2)-int(pos2)
				elif chr2==chr1_2 and terinfo[-1]=='5' and terinfo_2[0]=='5':
					i_dist2=int(pos2)-int(pos1_2)
				if i_dist2 < s_ra and i_dist2 > -1*SA_bp_v and i_dist2  < c_dist2:
					candidate_info2=[chr1_2,int(pos1_2),chr2_2,int(pos2_2),mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,str(allfrag_2),str(spfrag_2),str(safrag_2),ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2]
					c_dist2=i_dist2

				if chr2==chr2_2 and terinfo[-1]=='3' and terinfo_2[-1]=='3':
					i_dist2=int(pos2_2)-int(pos2)
				elif chr2==chr2_2 and terinfo[-1]=='5' and terinfo_2[-1]=='5':
					i_dist2=int(pos2)-int(pos2_2)
				if i_dist2 < s_ra and i_dist2 > -1*SA_bp_v and i_dist2  < c_dist2:
					candidate_info2=[chr1_2,int(pos1_2),chr2_2,int(pos2_2),mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,str(allfrag_2),str(spfrag_2),str(safrag_2),ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2]
					c_dist2=i_dist2
		if candidate_info1 != '' and candidate_info2 != '':
			output_list.append([chr1,int(pos1),chr2,int(pos2),mh,terinfo,svtype,precise,class1,'removal',str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
			output_id_list.append('\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype]))
			candidate_info1[9]='merged'
			new_allfrag=int(allfrag)-int(spfrag)
			candidate_info1[10]=str(int(candidate_info1[10])+int(new_allfrag))
			idx1='\t'.join([candidate_info1[0],str(candidate_info1[1]),candidate_info1[2],str(candidate_info1[3]),candidate_info1[4],candidate_info1[5],candidate_info1[6]])
			if idx1 not in output_id_list:
				output_list.append(candidate_info1)
				output_id_list.append(idx1)
			candidate_info2[9]='merged'
			candidate_info2[10]=str(int(candidate_info2[10])+int(new_allfrag))
			idx2='\t'.join([candidate_info2[0],str(candidate_info2[1]),candidate_info2[2],str(candidate_info2[3]),candidate_info2[4],candidate_info2[5],candidate_info2[6]])
			if idx2 not in output_id_list:
				output_list.append(candidate_info2)
				output_id_list.append(idx2)
		elif candidate_info1 != '' and unknown_bp2 == 'T':
			output_list.append([chr1,int(pos1),chr2,int(pos2),mh,terinfo,svtype,precise,class1,'removal',str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
			output_id_list.append('\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype]))
			candidate_info1[9]='merged'
			new_allfrag=int(allfrag)-int(spfrag)
			candidate_info1[10]=str(int(candidate_info1[10])+int(new_allfrag))
			idx1='\t'.join([candidate_info1[0],str(candidate_info1[1]),candidate_info1[2],str(candidate_info1[3]),candidate_info1[4],candidate_info1[5],candidate_info1[6]])
			if idx1 not in output_id_list:
				output_list.append(candidate_info1)
				output_id_list.append(idx1)
		elif unknown_bp1 == 'T' and candidate_info2 != '':
			output_list.append([chr1,int(pos1),chr2,int(pos2),mh,terinfo,svtype,precise,class1,'removal',str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
			output_id_list.append('\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype]))
			candidate_info2[9]='merged'
			new_allfrag=int(allfrag)-int(spfrag)
			candidate_info2[10]=str(int(candidate_info2[10])+int(new_allfrag))
			idx2='\t'.join([candidate_info2[0],str(candidate_info2[1]),candidate_info2[2],str(candidate_info2[3]),candidate_info2[4],candidate_info2[5],candidate_info2[6]])
			if idx2 not in output_id_list:
				output_list.append(candidate_info2)
				output_id_list.append(idx2)
		else:
			output_list.append([chr1,int(pos1),chr2,int(pos2),mh,terinfo,svtype,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
			output_id_list.append('\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype]))


# adding remained (not merged) precise SVs in initial_list to output_list.
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in initial_list:
	if precise == 'Y' and '\t'.join([chr1,pos1,chr2,pos2,mh,terinfo,svtype]) not in output_id_list:
		output_list.append([chr1,int(pos1),chr2,int(pos2),mh,terinfo,svtype,precise,class1,conc1,str(allfrag),str(spfrag),str(safrag),ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])


# Removal of both BP duplicates with high mh removal
output_2_list=[]
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_list:
	pos1=str(pos1);pos2=str(pos2)
	if mh == '.':
		output_2_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
	else:
		if class1 == 'initial':
			bp1_score='1000'; bp2_score='1000'
		else:
			bp1_score=int((class1.split('(')[1]).split(',')[0])
			bp2_score=int((class1.split(',')[1]).split(')')[0])
		bp1_dup="N";bp2_dup="N"
		bp1_score_list=[0]; bp2_score_list=[0]; mh1_list=[]; mh2_list=[];bp1_frag_list=[];bp2_frag_list=[]
		for [chr1_2,pos1_2,chr2_2,pos2_2,mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,allfrag_2,spfrag_2,safrag_2,ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2] in output_list:
			pos1_2=str(pos1_2);pos2_2=str(pos2_2)
			if mh_2 != '.':
				if chr1==chr1_2 and pos1 == pos1_2 and terinfo[0] == terinfo_2[0] and (chr2 != chr2_2 or pos2 != pos2_2 or terinfo[-1] != terinfo_2[-1]):
					bp1_dup = 'Y'
					mh1_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split(',')[1]).split(')')[0])
					bp2_score_list.append(score)
					bp2_frag_list.append(int(allfrag_2))
				elif chr1==chr2_2 and pos1 == pos2_2 and terinfo[0] == terinfo_2[-1] and (chr2 != chr1_2 or pos2 != pos1_2 or terinfo[-1] != terinfo_2[0]):
					bp1_dup = 'Y'
					mh1_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split('(')[1]).split(',')[0])
					bp2_score_list.append(score)
					bp2_frag_list.append(int(allfrag_2))
				elif chr2 == chr2_2 and pos2 == pos2_2 and terinfo[-1] == terinfo_2[-1] and (chr1 != chr1_2 or pos1 != pos1_2 or terinfo[0] != terinfo_2[0]):
					bp2_dup = 'Y'
					mh2_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split('(')[1]).split(',')[0])
					bp1_score_list.append(score)
					bp1_frag_list.append(int(allfrag_2))
				elif chr2 == chr1_2 and pos2 == pos1_2 and terinfo[-1] == terinfo_2[0] and (chr1 != chr2_2 or pos1 != pos2_2 or terinfo[0] != terinfo_2[-1]):
					bp2_dup = 'Y'
					mh2_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split(',')[1]).split(')')[0])
					bp1_score_list.append(score)
		if (bp1_dup == 'Y' and bp2_dup == 'Y' and abs(int(mh)) > min(mh1_list) and abs(int(mh))> min(mh2_list)):
			'blank'
		else:
			output_2_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])

# BP1 duplicate removal
output_3_list=[]
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_2_list:
	pos1=str(pos1);pos2=str(pos2)
	if mh == '.':
		output_3_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
	else:
		if class1 == 'initial':
			bp1_score='1000'; bp2_score='1000'
		else:
			bp1_score=int((class1.split('(')[1]).split(',')[0])
			bp2_score=int((class1.split(',')[1]).split(')')[0])
		bp1_dup="N";bp2_dup="N"
		bp1_score_list=[0]; bp2_score_list=[0]; mh1_list=[]; mh2_list=[];bp1_frag_list=[];bp2_frag_list=[]
		for [chr1_2,pos1_2,chr2_2,pos2_2,mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,allfrag_2,spfrag_2,safrag_2,ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2] in output_2_list:
			pos1_2=str(pos1_2);pos2_2=str(pos2_2)
			if mh_2 != '.':
				if chr1==chr1_2 and pos1 == pos1_2 and terinfo[0] == terinfo_2[0] and (chr2 != chr2_2 or pos2 != pos2_2 or terinfo[-1] != terinfo_2[-1]):
					bp1_dup = 'Y'
					mh1_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split(',')[1]).split(')')[0])
					bp2_score_list.append(score)
					bp2_frag_list.append(int(allfrag_2))
				elif chr1==chr2_2 and pos1 == pos2_2 and terinfo[0] == terinfo_2[-1] and (chr2 != chr1_2 or pos2 != pos1_2 or terinfo[-1] != terinfo_2[0]):
					bp1_dup = 'Y'
					mh1_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split('(')[1]).split(',')[0])
					bp2_score_list.append(score)
					bp2_frag_list.append(int(allfrag_2))
				elif chr2 == chr2_2 and pos2 == pos2_2 and terinfo[-1] == terinfo_2[-1] and (chr1 != chr1_2 or pos1 != pos1_2 or terinfo[0] != terinfo_2[0]):
					bp2_dup = 'Y'
					mh2_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split('(')[1]).split(',')[0])
					bp1_score_list.append(score)
					bp1_frag_list.append(int(allfrag_2))
				elif chr2 == chr1_2 and pos2 == pos1_2 and terinfo[-1] == terinfo_2[0] and (chr1 != chr2_2 or pos1 != pos2_2 or terinfo[0] != terinfo_2[-1]):
					bp2_dup = 'Y'
					mh2_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split(',')[1]).split(')')[0])
					bp1_score_list.append(score)

		if bp1_dup == 'Y' and bp2_dup == 'N' and bp2_score < max(bp2_score_list):
			'blank'
		elif bp1_dup == 'Y' and bp2_dup == 'N' and bp2_score == max(bp2_score_list):
			output_3_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, 'T'])
		else:
			output_3_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])

#BP2 duplicate removal
output_4_list=[]
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_3_list:
	pos1=str(pos1);pos2=str(pos2)
	if mh == '.':
		output_4_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
	else:
		if class1 == 'initial':
			bp1_score='1000'; bp2_score='1000'
		else:
			bp1_score=int((class1.split('(')[1]).split(',')[0])
			bp2_score=int((class1.split(',')[1]).split(')')[0])
		bp1_dup="N";bp2_dup="N"
		bp1_score_list=[0]; bp2_score_list=[0]; mh1_list=[]; mh2_list=[];bp1_frag_list=[];bp2_frag_list=[]
		for [chr1_2,pos1_2,chr2_2,pos2_2,mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,allfrag_2,spfrag_2,safrag_2,ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2] in output_3_list:
			pos1_2=str(pos1_2);pos2_2=str(pos2_2)
			if mh_2 != '.':
				if chr1==chr1_2 and pos1 == pos1_2 and terinfo[0] == terinfo_2[0] and (chr2 != chr2_2 or pos2 != pos2_2 or terinfo[-1] != terinfo_2[-1]):
					bp1_dup = 'Y'
					mh1_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split(',')[1]).split(')')[0])
					bp2_score_list.append(score)
					bp2_frag_list.append(int(allfrag_2))
				elif chr1==chr2_2 and pos1 == pos2_2 and terinfo[0] == terinfo_2[-1] and (chr2 != chr1_2 or pos2 != pos1_2 or terinfo[-1] != terinfo_2[0]):
					bp1_dup = 'Y'
					mh1_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split('(')[1]).split(',')[0])
					bp2_score_list.append(score)
					bp2_frag_list.append(int(allfrag_2))
				elif chr2 == chr2_2 and pos2 == pos2_2 and terinfo[-1] == terinfo_2[-1] and (chr1 != chr1_2 or pos1 != pos1_2 or terinfo[0] != terinfo_2[0]):
					bp2_dup = 'Y'
					mh2_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split('(')[1]).split(',')[0])
					bp1_score_list.append(score)
					bp1_frag_list.append(int(allfrag_2))
				elif chr2 == chr1_2 and pos2 == pos1_2 and terinfo[-1] == terinfo_2[0] and (chr1 != chr2_2 or pos1 != pos2_2 or terinfo[0] != terinfo_2[-1]):
					bp2_dup = 'Y'
					mh2_list.append(abs(int(mh_2)))
					if class1_2 == 'initial':
						score=1000
					elif class1_2 != 'initial':
						score=int((class1_2.split(',')[1]).split(')')[0])
					bp1_score_list.append(score)

		if bp1_dup == 'N' and bp2_dup == 'Y' and bp1_score < max(bp1_score_list):
			'blank'
		elif bp1_dup == 'N' and bp2_dup == 'Y' and bp1_score == max(bp1_score_list):
			output_4_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, 'T', unknown_bp2])
		else:
			output_4_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])

#filtering SAits without close other BP and newSAmate with initial non-precise SV
output_5_list=[]
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_4_list:
	if conc1 == 'removal' or (unknown_bp1 == 'T' and unknown_bp2 == 'T'):
		'blank'
	else:
		output_5_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
		
output_6_list=[]
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_5_list:
	if 'newSAmate' in class1:
		ter1=terinfo[0]
		ter2=terinfo[-1]
		initialN='F'
		for [chr1_2,pos1_2,chr2_2,pos2_2,mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,allfrag_2,spfrag_2,safrag_2,ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2] in output_5_list:
			if dellyid == dellyid_2 and class1_2 == 'initial':
				initialN='T'
				break
		if initialN =='F':
			output_6_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
	else:	
		output_6_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])

output_7_list=[]
close_ra = 120
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_6_list:
	if 'SAits' in class1:
		ter1=terinfo[0]
		ter2=terinfo[-1]
		bp1_close_oppo='F';bp2_close_oppo='F'
		for [chr1_2,pos1_2,chr2_2,pos2_2,mh_2,terinfo_2,svtype_2,precise_2,class1_2,conc1_2,allfrag_2,spfrag_2,safrag_2,ref1_2,ref2_2,mapq1_2,mapq2_2,tbginfo1_2,tbginfo2_2,nbginfo1_2,nbginfo2_2,pos_dist1_2,pos_dist2_2,dellyid_2, unknown_bp1_2, unknown_bp2_2] in output_6_list:
			if chr1==chr1_2 and pos1 == pos1_2 and chr2 == chr2_2 and pos2 == pos2_2 and terinfo == terinfo_2:  # pass same SV
				continue
			t1=terinfo_2[0]
			t2=terinfo_2[-1]
			if chr1==chr1_2 and ((ter1 == '3' and t1 == '5' and int(pos1)-int(pos1_2) > 0 and int(pos1)-int(pos1_2) < close_ra) or (ter1 == '5' and t1 == '3' and int(pos1_2)-int(pos1) > 0 and int(pos1_2) -int(pos1) < close_ra)):
				bp1_close_oppo='T'
			elif chr1==chr2_2 and ((ter1 == '3' and t2 == '5' and int(pos1)-int(pos2_2) > 0 and int(pos1)-int(pos2_2) < close_ra) or (ter1 == '5' and t2 == '3' and int(pos2_2)-int(pos1) > 0 and int(pos2_2) - int(pos1) < close_ra)):
				bp1_close_oppo='T'
			if chr2==chr1_2 and ((ter2 == '3' and t1 == '5' and int(pos2)-int(pos1_2) > 0 and int(pos2)-int(pos1_2) < close_ra) or (ter2 == '5' and t1 == '3' and int(pos1_2)-int(pos2) > 0 and int(pos1_2) -int(pos2) < close_ra)):
				bp2_close_oppo='T'
			elif chr2==chr2_2 and ((ter2 == '3' and t2 == '5' and int(pos2)-int(pos2_2) > 0 and int(pos2)-int(pos2_2) < close_ra) or (ter2 == '5' and t2 == '3' and int(pos2_2)-int(pos2) > 0 and int(pos2_2) - int(pos2) < close_ra)):
				bp2_close_oppo='T'
		if bp1_close_oppo=='T' or bp2_close_oppo=='T':
			output_7_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
	else:
		output_7_list.append([chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2])
		
# filtering removal-marked SV and editing and deduplicates of unknown BP.
output_8_list=[]
output_8_id_list=[]
for [chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2] in output_7_list:
	chr1n=change_chr_to_int(chr1)
	chr2n=change_chr_to_int(chr2)
	if unknown_bp1 == 'F' and unknown_bp2 == 'F':
		idx=[chr1n, int(pos1), chr2n, int(pos2),chr1,pos1,chr2,pos2,mh,terinfo,svtype]
		info_list=[chr1n, int(pos1), chr2n, int(pos2),chr1,pos1,chr2,pos2,mh,terinfo,svtype,precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
	elif unknown_bp1=='T':
		idx=[chr2n, int(pos2), 0,0, chr2,pos2,'.','.','.',terinfo[-1]+'to.','.']
		info_list=[chr2n, int(pos2), 0,0,chr2,pos2,'.','.','.',terinfo[-1]+'to.','.',precise,class1,conc1,allfrag,spfrag,safrag,ref2,ref1,mapq2,mapq1,tbginfo2, tbginfo1, nbginfo2, nbginfo1, pos_dist2, pos_dist1, dellyid, unknown_bp2, unknown_bp1]
	elif unknown_bp2 == 'T':
		idx=[chr1n, int(pos1),0,0, chr1,pos1,'.','.','.',terinfo[0]+'to.','.']
		info_list=[chr1n, int(pos1),0,0,chr1,pos1,'.','.','.',terinfo[0]+'to.','.',precise,class1,conc1,allfrag,spfrag,safrag,ref1, ref2, mapq1, mapq2, tbginfo1, tbginfo2, nbginfo1, nbginfo2, pos_dist1, pos_dist2, dellyid, unknown_bp1, unknown_bp2]
	if idx not in output_8_id_list:
		output_8_list.append(info_list)
		output_8_id_list.append(idx)
print(len(output_8_list))

output_8_list.sort(key=itemgetter(0,1,2,3))
for info in output_8_list:
	out_file.write('\t'.join(info[4:])+'\n')

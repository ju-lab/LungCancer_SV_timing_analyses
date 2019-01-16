
import sys
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.fi','w')
n=0;m=0
in_line=in_file.readline().strip()
while in_line:
        if in_line[0:4]=='#CHR':
                out_file.write(in_line+'\n')
	elif in_line[0]=='#':
		'blank'
        else:
		n+=1
                in_indi=in_line.split('\t')
		#column number assign
		ninfo=in_indi[10]
                npinfo=in_indi[11]
                chr1=in_indi[14]
		pos1=in_indi[15]
		chr2=in_indi[16]
		pos2=in_indi[17]
		terinfo=in_indi[19]
		svtype=in_indi[20]
		if svtype=='INS':
			in_line=in_file.readline().strip()
			continue
#		mapq=in_indi[21]
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
		pos1=int(pos1);pos2=int(pos2)
		dist=abs(pos2-pos1)
		nprc=int(npinfo.split(';')[0])
		npsn=int(npinfo.split(';')[1])
		ngt=ninfo.split(':')[0]
		ter1=terinfo.split('to')[0]
		ter2=terinfo.split('to')[1]

		allfrag=int(tfinfo.split(';')[2])
		spfrag=int(tfinfo.split(';')[3])
		safrag=int(tfinfo.split(';')[4])
		vaf1=(tfinfo.split(';')[5][:-1])
		if vaf1=='N':
			vaf1=0
		else:
			vaf1=float(vaf1)
		vaf2=(tfinfo.split(';')[6][:-1])
		if vaf2=='N':
			vaf2=0
		else:
			vaf2=float(vaf2)
		pnfrag=int(nfinfo.split(';')[0])
		pnsa=int(nfinfo.split(';')[1])
		pnfc1=int(nfinfo.split(';')[2])
		pnfc2=int(nfinfo.split(';')[3])

		tbgcn1=int(tbginfo1.split(';')[3])
		tbgbin1=int(tbginfo1.split(';')[4])
		tbgcnt1=int(tbginfo1.split(';')[5])
		tbgbint1=int(tbginfo1.split(';')[6])
		tbgcn2=int(tbginfo2.split(';')[3])
		tbgbin2=int(tbginfo2.split(';')[4])
		tbgcnt2=int(tbginfo2.split(';')[5])
		tbgbint2=int(tbginfo2.split(';')[6])
		
		nbgcn1=int(nbginfo1.split(';')[3])
		nbgbin1=int(nbginfo1.split(';')[4])
		nbgcnt1=int(nbginfo1.split(';')[5])
		nbgbint1=int(nbginfo1.split(';')[6])
		nbgcn2=int(nbginfo2.split(';')[3])
		nbgbin2=int(nbginfo2.split(';')[4])
		nbgcnt2=int(nbginfo2.split(';')[5])
		nbgbint2=int(nbginfo2.split(';')[6])
	
		pnsc=int(pnsc)
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

		if pos_dist1=='NA':
			max_pos1=1000
		else:
			max_pos1=int(pos_dist1.split(';')[2])
		if pos_dist2=='NA':
			max_pos2=1000
		else:
			max_pos2=int(pos_dist2.split(';')[2])
	
		new1sa_max=0;sa_mate_n1=0;tsan1=0;tsan2=0   # tsan = target SA number
		if new1 != 'NA':
			target='off'
			tom_list=[] # target other mate list
			osan_list=[]  # other SA number list
			otsan=0 
			for info in new1.split(','):
				san=int((info.split('(')[1]).split(')')[0])
				c1=info.split(':')[0]
				p1=(info.split(';')[0]).split(':')[1]
				c2=(info.split(';')[1]).split(':')[0]
				p2=(info.split(';')[1]).split(':')[1]
				t1=(info.split(';')[4]).split('to')[0]
				t2=(info.split(';')[4]).split('to')[1]
				p1=int(p1); p2=int(p2)
				if c1==chr1 and t1==ter1 and abs(pos1-p1) < 5:
					if c2 == chr2 and t2==ter2 and abs(pos2-p2) <5:
						tsan1+=san
						target='on'
					else:
						tom_list.append(';'.join([c2,str(p2),t2]))
				elif c2==chr1 and t2==ter1 and abs(pos1-p2) < 5:
					if c1 == chr2 and t1 == ter2 and abs(p1-pos2) < 5:
						tsan1 += san
						target='on'
					else:
						tom_list.append(';'.join([c1,str(p1),t1]))
				if target == 'off':
					osan_list.append(san)
			if len(osan_list)==0:
				new1sa_max=0
			else:
				new1sa_max=max(osan_list)
		new2sa_max=0;sa_mate_n2=0
		if new2 != 'NA':
			target='off'
			tom_list=[]
			osan_list=[]
			for info in new2.split(','):
				san=int((info.split('(')[1]).split(')')[0])
				c1=info.split(':')[0]
				p1=(info.split(';')[0]).split(':')[1]
				c2=(info.split(';')[1]).split(':')[0]
				p2=(info.split(';')[1]).split(':')[1]
				t1=(info.split(';')[4]).split('to')[0]
				t2=(info.split(';')[4]).split('to')[1]
				p1=int(p1); p2=int(p2)
				if c1==chr2 and t1==ter2 and abs(pos2-p1) < 5:
					if c2 == chr1 and t2 == ter1 and abs(p2-pos1) < 5:
						tsan2 += san
						target='on'
					else:
						tom_list.append(';'.join([c2,str(p2),t2]))
				elif c2==chr2 and t2==ter2 and abs(pos2-p2) < 5:
					if c1 == chr1 and t1 == ter1 and abs(p1-pos1) < 5:
						tsan2 += san
						target='on'
					else:
						tom_list.append(';'.join([c1,str(p1),t1]))
				if target == 'off':
					osan_list.append(san)
			if len(osan_list)==0:
				new2sa_max=0
			else:
				new2sa_max=max(osan_list)

		# breakpoint probelm definition
		bp1_problem= 'F'; bp2_problem = 'F'
		if med_mapq1 < 40 or (max_pos1 < 2 and new1sa_max < 2) or (dist >=1500 and nbgbin1 >=5):
			bp1_problem = 'T'
		if med_mapq2 < 40 or (max_pos2 < 2 and new2sa_max < 2) or (dist >=1500 and nbgbin2 >=5):
			bp2_problem = 'T'
	
		# breakpoint failure definition and fragment number readjustment.
		bp1_failure='F'; bp2_failure = 'F'
		if med_mapq1 <5 or nbgbin1 >= 50:
			bp1_failure = 'T'
		if med_mapq2 <5 or nbgbin2 >= 50:
			bp2_failure = 'T'
		if (bp1_failure =='F' and bp2_failure == 'T' and tsan1 == 0 ) or (bp1_failure == 'T' and bp2_failure == 'F' and tsan2 ==0):
			allfrag = allfrag - spfrag
			spfrag=0
			safrag=0
			
		#filter criteria
		if pnsa > 0 :
			'blank'
		elif (svtype != 'TRA' and dist < 1000 and npsn >=2) or ((svtype == 'TRA' or dist >=1000) and npsn >= 1):
			'blank'
		elif ((svtype != 'DEL' and svtype != 'INV') or dist >= 1000) and pnfrag > 0:
			'blank'
		elif med_mapq1 < 40 and med_mapq2 < 40 or min(max_mapq1, max_mapq2) < 20:
			'blank'
		elif min(nbgbin1, nbgbin2) >= 5 and (svtype == 'TRA' or dist >= 1500):
			'blank'
		elif max(new1sa_max, new2sa_max) < 2 and max(max_pos1, max_pos2) < 2:
			'blank'
		elif bp1_problem == 'T' and bp2_problem == 'T':
			'blank'
		elif allfrag < 10 and safrag <2 and min(nbgbin1, nbgbin2) >=3 and max(nbgbin1, nbgbin2) >=10 and (svtype == 'TRA' or dist > 1500):
			'blank'
		elif allfrag < 10 and safrag == 0 and ((nbgbin1 >=3 and med_mapq2 < 30) or (nbgbin2 >=3 and med_mapq1 < 30)):
			'blank'
		elif allfrag < 5 and safrag == 0 and (min(med_mapq1, med_mapq2) < 40 or (max_pos1 < 2 and new1sa_max <2) or (max_pos2 <2 and new2sa_max <2)):
			'blank'
		elif svtype == 'DUP' and dist < 10000 and safrag == 0:
			'blank'
		elif svtype == 'INV' and dist < 5000 and allfrag <5:
			'blank'
		elif svtype != 'TRA' and dist < 1000 and allfrag <5:
			'blank'
		elif (svtype == 'DEL' or svtype == 'INV') and dist < 1000 and pnfrag > 0 and spfrag <= 1 :
			'blank'
		elif pnsc > 0 and allfrag-spfrag < 3 and nbgbint1 < 10 and nbgbint2 < 10:
			'blank'
		elif vaf1 < 3 and vaf2 < 3: 
			'blank'
		elif svtype == 'DEL' and dist < 50:  #optional!!
			'blank'
		elif allfrag < 3:  # optional!!
			'blank'
		else:
			'blank'
			m+=1
			out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()
print(n)
print(m)

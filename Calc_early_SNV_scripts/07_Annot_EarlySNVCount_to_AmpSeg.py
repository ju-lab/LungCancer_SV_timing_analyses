#Arg1: CNV file (edit.boradAMP)
#Arg2: SNV file (.proba)
#Arg3: sampleid
#Arg4: outdir

import sys, math
print('###'+sys.argv[0])
print(sys.argv[1])
seqz_file=open(sys.argv[1])  
snv_file=open(sys.argv[2])    
out_file=open(sys.argv[4]+'/'+sys.argv[3]+'.ampseg.earlysnv.txt','w')  #id
header_list=['#chr', 'pos1', 'pos2', 'total_cn','major_cn', 'minor_cn', 'broadAMP','dist','early_SNV','early_SNV_range','early_SNVrate/Mb', 'early_SNVrate_CI','early_SNV_var']
out_file.write('\t'.join(header_list)+'\n')
def count_earlySNV(chr1, pos1, pos2):
	pos1=int(pos1); pos2=int(pos2)
	snv_file.seek(0)
	snv_line=snv_file.readline().strip()
	early=0;earlyV=0
	while snv_line:
		if snv_line[0]=='#':
			'blank'
		else:
			snv_indi=snv_line.split('\t')
			snv_chr=snv_indi[0]
			snv_pos=int(snv_indi[1])
			kat=snv_indi[16]
			pclonal=snv_indi[17]
			pearly=snv_indi[18]
			plate=snv_indi[19]
			psubcl=snv_indi[20]
			if snv_chr == chr1 and snv_pos >=pos1 and snv_pos <=pos2 and kat=='F':
				if (float(pclonal)==0 or math.isnan(float(pclonal))==True) and pearly=='.':
					pearly=0
				elif float(pclonal)==1 and pearly=='.':
					pearly=1
				try:
					early += float(pearly)
				except:
					print(snv_line)
				earlyV += float(pearly)*(1-float(pearly))
		snv_line=snv_file.readline().strip()
	learly=round(early-(1.96*math.sqrt(earlyV)),4)
	hearly=round(early+(1.96*math.sqrt(earlyV)),4)
	early = round(early,4)
	earlyV= round(earlyV,4)
	return [early, learly, hearly, earlyV]

snv_file.seek(0)
snv_line=snv_file.readline().strip()
tcln=0;tclnV=0
tsnvn=0
while snv_line:
	if snv_line[0]=='#':
		'blank'
	else:
		snv_indi=snv_line.split('\t')
		tsnvn +=1
		if snv_indi[17]=='.':
			'blank'
		else:
			pclonal=float(snv_indi[17])
			kat=snv_indi[16]
			if kat == 'F':
				tcln += pclonal
				tclnV += pclonal*(1-pclonal)
	snv_line=snv_file.readline().strip()
ltcln=round(tcln-(1.96*math.sqrt(tclnV)),4)
htcln=round(tcln+(1.96*math.sqrt(tclnV)),4)
tcln=round(tcln,4)
tclr=round(tcln*1000000/float(2800000000),4)
ltclr=round(ltcln*1000000/float(2800000000),4)
htclr=round(htcln*1000000/float(2800000000),4)
tsnvr=round(tsnvn*1000000/float(2800000000),4)

seqz_line=seqz_file.readline().strip()
while seqz_line:
	if seqz_line[0]=='#':
		'blank'
	else:
		seqz_indi=seqz_line.split('\t')
		chr1=seqz_indi[0]
		pos1=int(seqz_indi[1])
		pos2=int(seqz_indi[2])
		dist=pos2-pos1
		totCN=int(seqz_indi[3])
		MCN=int(seqz_indi[4])
		mCN=int(seqz_indi[5])
		if MCN >= 2:
			ec_list=count_earlySNV(chr1,pos1,pos2)
			earlyn=ec_list[0]
			learlyn=ec_list[1]
			hearlyn=ec_list[2]
			earlyvar=ec_list[3]
			if mCN <= MCN*0.75:
				earlyr=round(2*earlyn*1000000/float(dist),4)
				learlyr=round(2*learlyn*1000000/float(dist),4)
				hearlyr=round(2*hearlyn*1000000/float(dist),4)
			elif mCN > MCN*0.75:
				earlyr=round(earlyn*1000000/float(dist),4)
				learlyr=round(learlyn*1000000/float(dist),4)
				hearlyr=round(hearlyn*1000000/float(dist),4)
			if tclr==0:
				earlyar='NA'; learlyar='NA';hearlyar='NA'
			else:
				earlyar=round(earlyr/float(tclr),4)
				learlyar=round(learlyr/float(tclr),4)
				hearlyar=round(hearlyr/float(tclr),4)
			if tsnvr == 0:
				earlyar2='NA'; learlyar2='NA'; hearlyar2='NA'
			else:
				earlyar2=round(earlyr/float(tsnvr),4)
				learlyar2=round(learlyr/float(tsnvr),4)
				hearlyar2=round(hearlyr/float(tsnvr),4)
			info_list=seqz_indi[0:7]+[str(dist),str(earlyn),str(learlyn)+':'+str(hearlyn),str(earlyr),str(learlyr)+':'+str(hearlyr),str(earlyvar) ]
			out_file.write('\t'.join(info_list)+'\n')
		else:
			info_list=seqz_indi[0:7]+[str(dist),'.','.','.','.','.']
			out_file.write('\t'.join(info_list)+'\n')
	seqz_line=seqz_file.readline().strip()

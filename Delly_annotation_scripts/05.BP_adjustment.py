import sys
print('### BP adjustment')
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.BPadj','w')

tcol=13  # column number of tumor BP info (start from 1)
ncol=14  # normal BP info(start from 1)
delly_tcol=10 # delly tumor info(start from 1)

tcol -=1
ncol -=1
delly_tcol -=1

in_line=in_file.readline().strip()
while in_line:
	if in_line[0:4]=='#CHR':
		out_file.write(in_line+'\tre_chr1\tre_pos1\tre_chr2\tre_pos2\tMH\tterminal\tSVtype\tMAPQ\tDV\tRV\ttSA\tnSA\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		tBP=in_indi[tcol]
		nBP=in_indi[ncol]
		chr1=in_indi[0]
		pos1=int(in_indi[1])
		chr2=(in_indi[7].split('CHR2=')[1]).split(';')[0]
		pos2=int((in_indi[7].split('END=')[1]).split(';')[0])
		terinfo=(in_indi[7].split('CT=')[1]).split(';')[0]
		mapq=(in_indi[7].split('MAPQ=')[1]).split(';')[0]
		svtype=in_indi[2][0:3]
		DV=str(int(in_indi[delly_tcol].split(':')[-3]))
		RV=str(int(in_indi[delly_tcol].split(':')[-1]))
		MH='.'
		tSA=0;nSA=0
		if tBP=='NA':
			'blank'
		else:
			tBP_indi=tBP.split(',')
			current_tSA=0
			candidate_list=[]
			for BP in tBP_indi:
				if terinfo in BP and svtype in BP:
					tSA_count=int(BP.split('(')[1][:-1])
					if tSA_count > current_tSA:
						candidate_list=[]
						candidate_list.append(BP)
						current_tSA=tSA_count
					elif tSA_count == current_tSA:
						candidate_list.append(BP)
					elif tSA_count < current_tSA:
						'blank'
			if len(candidate_list)==1:
				final_tBP=candidate_list[0]
				final_tBP=final_tBP.replace('NC_','NC')
				bpchr1=(final_tBP.split('_')[0]).split(':')[0]
				bppos1=int((final_tBP.split('_')[0]).split(':')[1])
				bpchr2=(final_tBP.split('_')[1]).split(':')[0]
				bppos2=int((final_tBP.split('_')[1]).split(':')[1])
				MH=str(int(final_tBP.split('_')[2]))
				tSA=current_tSA
				bpchr1=bpchr1.replace('NC','NC_'); bpchr2=bpchr2.replace('NC','NC_')

				if chr1!=bpchr1:
					print('chromosome1 error')
					print(candidate_list)
					sys.exit(1)
				if chr2!=bpchr2:
					print('chromosome2 error')
					print(in_line)
					sys.exit(1)

				pos1=bppos1
				pos2=bppos2
			elif len(candidate_list) >1:
				current_distsum=99999999999
				for BP in candidate_list:
					BP=BP.replace('NC_','NC')
					bpchr1=(BP.split('_')[0]).split(':')[0]
					bppos1=int((BP.split('_')[0]).split(':')[1])
					bpchr2=(BP.split('_')[1]).split(':')[0]
					bppos2=int((BP.split('_')[1]).split(':')[1])
					bpMH=str(int(BP.split('_')[2]))
					
					bpchr1=bpchr1.replace('NC','NC_'); bpchr2=bpchr2.replace('NC','NC_')
					if chr1 != bpchr1:
						print('chromosome1 error')
						print(in_line)
						sys.exit(1)
					if chr2 != bpchr2:
						print('chromosome2 error')
						print(in_line)
						sys.exit(1)
					dist1=abs(pos1-bppos1)
					dist2=abs(pos2-bppos2)
					distsum=dist1+dist2
					if distsum < current_distsum:
						final_info=[bpchr1,bppos1,bpchr2,bppos2,bpMH, current_tSA]
						current_distsum=distsum
				chr1=final_info[0];pos1=final_info[1];chr2=final_info[2];pos2=final_info[3];MH=final_info[4];tSA=final_info[5]
		if nBP == 'NA':
			'blank'
		else:
			nBP_indi=nBP.split(',')	
			nSA_candidate=[]
			for BP in nBP_indi:
				if terinfo in BP and svtype in BP:
					BP=BP.replace('NC_','NC')
					bpchr1=BP.split('_')[0]
					bppos1=int((BP.split('_')[0]).split(':')[1])
					bpchr2=(BP.split('_')[1]).split(':')[0]
					bppos2=int((BP.split('_')[1]).split(':')[1])
					bpMH=str(int(BP.split('_')[2]))
					nSA_count=int(BP.split('(')[1][:-1])

					bpchr1=bpchr1.replace('NC','NC_'); bpchr2=bpchr2.replace('NC','NC_')

					dist1=abs(pos1-bppos1)
					dist2=abs(pos2-bppos2)
					if dist1 <= 1 and dist2 <=1:
						nSA_candidate.append(nSA_count)
			if len(nSA_candidate)==0:
				nSA=0
			else:
				nSA=max(nSA_candidate)
		info_list=[chr1,str(pos1),chr2,str(pos2),MH,terinfo, svtype, mapq, DV, RV, str(tSA), str(nSA)]
		out_file.write(in_line+'\t'+'\t'.join(info_list)+'\n')
	in_line=in_file.readline().strip()

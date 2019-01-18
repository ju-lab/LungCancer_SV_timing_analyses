import sys,re
print(sys.argv[1])
in_file=open(sys.argv[1])  # *.annot
sampleid=sys.argv[2]
mttype=sys.argv[3]  # snv or indel
out_file=open(sampleid+'.'+mttype+'.fsa','w')
idx_list=[]
in_line=in_file.readline().strip()
n=0
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		if mttype == 'snv':
			annot=in_indi[57]
		elif mttype=='indel':
			annot=in_indi[59]
		else:
			print('ERROR: incorrect mttype')
			sys.exit(1)
		if 'nsSNP' in annot or 'frameshift' in annot or 'inframe' in annot:
			n+=1
		annot_indi=annot.split(';')
		for info in annot_indi:
			info_indi=info.split(',')
			if info_indi[0]=='nsSNP':
				altaa=(info_indi[6].split('>')[1]).split(']')[0]
				if altaa=='*':  #nonsense
					continue
				idx=info_indi[1]+'_'+info_indi[5]
				if idx not in idx_list:
					aa1=info_indi[6].split('[')[0].upper()
					if len(aa1)==9:
						aa1=aa1[1:]
					aa2=(info_indi[6].split('>')[1]).split(']')[0]
					aa3=info_indi[6].split(']')[1].upper().replace('*','')
					if len(aa3)==9:
						aa3=aa3[:-1]
					aaseq=aa1+aa2+aa3
					out_file.write('>s'+str(n)+'_'+idx+'\n')
					out_file.write(aaseq+'\n')
					idx_list.append(idx)
			elif info_indi[0]=='frameshift_CDS_INDEL':
				altaa=(info_indi[6].split('>')[1]).split(']')[0]
				if altaa=='*':	#nonsense
					continue
				sobj=re.search(r'([0-9]+)([A-Z])',info_indi[5])
				idx=info_indi[1]+'_'+sobj.group(2)+sobj.group(1)+'fs'
				if idx not in idx_list:
					aa1=info_indi[6].split('[')[0].upper()
					if len(aa1)==9:
						aa1=aa1[1:]
					aa2=(info_indi[6].split('>')[1]).split(']')[0]
					aa3=info_indi[6].split(']')[1].upper().replace('*','')
					aaseq=aa1+aa2+aa3
					out_file.write('>i'+str(n)+'_'+idx+'\n')
					out_file.write(aaseq+'\n')
					idx_list.append(idx)
			elif info_indi[0]=='inframe_CDS_INDEL':
				altaa=(info_indi[6].split('>')[1]).split(']')[0]
				if altaa=='*':	#nonsense
					continue
				if info_indi[5][-1]=='>':
					idx=info_indi[1]+'_'+info_indi[5].replace('>','del')
				else:
					idx=info_indi[1]+'_'+info_indi[5].replace('>','ins')
				if idx not in idx_list:
					aa1=info_indi[6].split('[')[0].upper()
					if len(aa1)==9:
						aa1=aa1[1:]
					aa2=(info_indi[6].split('>')[1]).split(']')[0]
					aa3=info_indi[6].split(']')[1].upper().replace('*','')
					if len(aa3)==9:
						aa3=aa3[:-1]
					aaseq=aa1+aa2+aa3
					out_file.write('>i'+str(n)+'_'+idx+'\n')
					out_file.write(aaseq+'\n')
					idx_list.append(idx)
	in_line=in_file.readline().strip()


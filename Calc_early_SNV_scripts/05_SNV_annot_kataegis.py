#This script was made 2018-06-14
import sys
print('###'+sys.argv[0])
print(sys.argv[1])
in_file=open(sys.argv[1]) # sorted snv
out_file=open(sys.argv[1]+'.kat','w')
in_line=in_file.readline().strip()
dist=1000   # kataegis distance criteria
pprv_chr=''
pprv_pos=''
prv_chr=''
prv_pos=''
while in_line:
	if in_line[0:6]=='#CHROM':
		out_file.write(in_line+'\tkataegis\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		kat='off'
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=int(in_indi[1])

		markpos=in_file.tell()
		next_line=in_file.readline().strip()
		if next_line=='': 
			next_chr=''; nnext_chr=''
		else:
			next_indi=next_line.split('\t')
			next_chr=next_indi[0]
			next_pos=int(next_indi[1])
			nnext_line=in_file.readline().strip()
			if nnext_line=='':
				nnext_chr=''
			else:
				nnext_indi=nnext_line.split('\t')
				nnext_chr=nnext_indi[0]
				nnext_pos=int(nnext_indi[1])
		in_file.seek(markpos)
		

		if pprv_chr != '':
			if chr1==pprv_chr and pos1-pprv_pos < dist:
				kat='on'
		if prv_chr != '' and next_chr !='':		
			if next_chr==prv_chr and next_pos-prv_pos < dist:
				kat='on'
		if nnext_chr != '':
			if nnext_chr==chr1 and nnext_pos-pos1 < dist:
				kat='on'
		
		if kat =='on':
			out_file.write(in_line+'\tT\n')
		elif kat == 'off':
			out_file.write(in_line+'\tF\n')

		pprv_chr=prv_chr; pprv_pos=prv_pos
		prv_chr=chr1; prv_pos=pos1

	in_line=in_file.readline().strip()

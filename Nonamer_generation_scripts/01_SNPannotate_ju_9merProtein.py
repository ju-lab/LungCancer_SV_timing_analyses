import sys

def blosum90(wtAA,varAA):
	wtAA=int(wtAA.upper().replace("A","1").replace("R","2").replace("N","3").replace("D","4").replace("C","5").replace("Q","6").replace("E","7").replace("G","8").replace("H","9").replace("I","10").replace("L","11").replace("K","12").replace("M","13").replace("F","14").replace("P","15").replace("S","16").replace("T","17").replace("W","18").replace("Y","19").replace("V","20").replace("B","21").replace("J","22").replace("Z","23").replace("X","24").replace("*","25"))
	varAA=int(varAA.upper().replace("A","1").replace("R","2").replace("N","3").replace("D","4").replace("C","5").replace("Q","6").replace("E","7").replace("G","8").replace("H","9").replace("I","10").replace("L","11").replace("K","12").replace("M","13").replace("F","14").replace("P","15").replace("S","16").replace("T","17").replace("W","18").replace("Y","19").replace("V","20").replace("B","21").replace("J","22").replace("Z","23").replace("X","24").replace("*","25"))
	
	blosum90matrix=[\
["0"," A"," R"," N"," D"," C"," Q"," E"," G"," H"," I"," L"," K"," M"," F"," P"," S"," T"," W"," Y"," V"," B"," J"," Z"," X"," *"],\
["A"," 5","-2","-2","-3","-1","-1","-1"," 0","-2","-2","-2","-1","-2","-3","-1"," 1"," 0","-4","-3","-1","-2","-2","-1","-1","-6"],\
["R","-2"," 6","-1","-3","-5"," 1","-1","-3"," 0","-4","-3"," 2","-2","-4","-3","-1","-2","-4","-3","-3","-2","-3"," 0","-1","-6"],\
["N","-2","-1"," 7"," 1","-4"," 0","-1","-1"," 0","-4","-4"," 0","-3","-4","-3"," 0"," 0","-5","-3","-4"," 5","-4","-1","-1","-6"],\
["D","-3","-3"," 1"," 7","-5","-1"," 1","-2","-2","-5","-5","-1","-4","-5","-3","-1","-2","-6","-4","-5"," 5","-5"," 1","-1","-6"],\
["C","-1","-5","-4","-5"," 9","-4","-6","-4","-5","-2","-2","-4","-2","-3","-4","-2","-2","-4","-4","-2","-4","-2","-5","-1","-6"],\
["Q","-1"," 1"," 0","-1","-4"," 7"," 2","-3"," 1","-4","-3"," 1"," 0","-4","-2","-1","-1","-3","-3","-3","-1","-3"," 5","-1","-6"],\
["E","-1","-1","-1"," 1","-6"," 2"," 6","-3","-1","-4","-4"," 0","-3","-5","-2","-1","-1","-5","-4","-3"," 1","-4"," 5","-1","-6"],\
["G"," 0","-3","-1","-2","-4","-3","-3"," 6","-3","-5","-5","-2","-4","-5","-3","-1","-3","-4","-5","-5","-2","-5","-3","-1","-6"],\
["H","-2"," 0"," 0","-2","-5"," 1","-1","-3"," 8","-4","-4","-1","-3","-2","-3","-2","-2","-3"," 1","-4","-1","-4"," 0","-1","-6"],\
["I","-2","-4","-4","-5","-2","-4","-4","-5","-4"," 5"," 1","-4"," 1","-1","-4","-3","-1","-4","-2"," 3","-5"," 3","-4","-1","-6"],\
["L","-2","-3","-4","-5","-2","-3","-4","-5","-4"," 1"," 5","-3"," 2"," 0","-4","-3","-2","-3","-2"," 0","-5"," 4","-4","-1","-6"],\
["K","-1"," 2"," 0","-1","-4"," 1"," 0","-2","-1","-4","-3"," 6","-2","-4","-2","-1","-1","-5","-3","-3","-1","-3"," 1","-1","-6"],\
["M","-2","-2","-3","-4","-2"," 0","-3","-4","-3"," 1"," 2","-2"," 7","-1","-3","-2","-1","-2","-2"," 0","-4"," 2","-2","-1","-6"],\
["F","-3","-4","-4","-5","-3","-4","-5","-5","-2","-1"," 0","-4","-1"," 7","-4","-3","-3"," 0"," 3","-2","-4"," 0","-4","-1","-6"],\
["P","-1","-3","-3","-3","-4","-2","-2","-3","-3","-4","-4","-2","-3","-4"," 8","-2","-2","-5","-4","-3","-3","-4","-2","-1","-6"],\
["S"," 1","-1"," 0","-1","-2","-1","-1","-1","-2","-3","-3","-1","-2","-3","-2"," 5"," 1","-4","-3","-2"," 0","-3","-1","-1","-6"],\
["T"," 0","-2"," 0","-2","-2","-1","-1","-3","-2","-1","-2","-1","-1","-3","-2"," 1"," 6","-4","-2","-1","-1","-2","-1","-1","-6"],\
["W","-4","-4","-5","-6","-4","-3","-5","-4","-3","-4","-3","-5","-2"," 0","-5","-4","-4","11"," 2","-3","-6","-3","-4","-1","-6"],\
["Y","-3","-3","-3","-4","-4","-3","-4","-5"," 1","-2","-2","-3","-2"," 3","-4","-3","-2"," 2"," 8","-3","-4","-2","-3","-1","-6"],\
["V","-1","-3","-4","-5","-2","-3","-3","-5","-4"," 3"," 0","-3"," 0","-2","-3","-2","-1","-3","-3"," 5","-4"," 1","-3","-1","-6"],\
["B","-2","-2"," 5"," 5","-4","-1"," 1","-2","-1","-5","-5","-1","-4","-4","-3"," 0","-1","-6","-4","-4"," 5","-5"," 0","-1","-6"],\
["J","-2","-3","-4","-5","-2","-3","-4","-5","-4"," 3"," 4","-3"," 2"," 0","-4","-3","-2","-3","-2"," 1","-5"," 4","-4","-1","-6"],\
["Z","-1"," 0","-1"," 1","-5"," 5"," 5","-3"," 0","-4","-4"," 1","-2","-4","-2","-1","-1","-4","-3","-3"," 0","-4"," 5","-1","-6"],\
["X","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-6"],\
["*","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6","-6"," 1"],\
]
	return int(blosum90matrix[wtAA][varAA])

def tRNAlocus(gidx,gstr,mutpos): #  tRNAinfo=tRNAlocus(gene_index,gene_strand,pos) ####FIXING
	tRNAinfofn="/home/users/ysj/Database/02_genes/01_UCSC_refgene/refGene_chrMTs_tRNAstructure.txt"
	tRNAfile=file(tRNAinfofn)
	tRNAline=tRNAfile.readline() #header
	tRNAline=tRNAfile.readline()

	while tRNAline:
		tRNAline_split=tRNAline.rstrip().split("\t")
		tRNAindex=tRNAline_split[1]
		if tRNAindex != gidx:
			tRNAline=tRNAfile.readline()
			continue
		break
	if tRNAline=="":
		print("Problem! No tRNA information: MT:"+str(mutpos)+"\t"+gidx)
		raw_input()
	tRNAstart=int(tRNAline_split[4]) # 0 base
	tRNAstop=int(tRNAline_split[5])
	tRNAstrand=tRNAline_split[3]
	tRNAstructure=tRNAline_split[-15:]
	tRNAannotation=["Acc-stem1","Acc-D_link","D-stem1","D-loop","D-stem2","D-Ac_link","Ac-stem1","Ac-loop","Ac-stem2","V-region","T-stem1","T-loop","T-stem2","Acc-stem2","TerminalBase"]

	if tRNAstrand!=gstr:
		print("Problem! tRNA strand is different between RefGene and tRNA locus file. "+gidx)
		raw_input()

	if tRNAstrand=="+":
		baseoffset=mutpos-tRNAstart
	else: # tRNAstrand=="-":
		baseoffset=tRNAstop-mutpos+1

	tRNAcursor=0
	thisAnnotation=""
	thisAnnotationContext=""
	for n, tstr in enumerate(tRNAstructure):
		tstr=int(tstr)
		if baseoffset <= (tRNAcursor + tstr) :
			thisAnnotation=tRNAannotation[n]
			if thisAnnotation!="Ac-loop":
				thisAnnotationContext= str(baseoffset-tRNAcursor) +"/"+str(tstr)
			else:
				thisAnnotationContext= (baseoffset-tRNAcursor)
				if thisAnnotationContext==3:
					thisAnnotationContext="3rdAntiCodon"
				elif thisAnnotationContext==4:
					thisAnnotationContext="2ndAntiCodon"
				elif thisAnnotationContext==5:
					thisAnnotationContext="1stAntiCodon"
				else:
					thisAnnotationContext= str(baseoffset-tRNAcursor) +"/"+str(tstr)
			break
		tRNAcursor+=tstr
	return [thisAnnotation, thisAnnotationContext]

def get_sequence(gs_chr,gs_start,gs_stop):
	seq_fn="/home/users/ysj/Database/01_reference/01_hg19/chr"+gs_chr+".fa"
	seq_file=file(seq_fn)
	seq_line=seq_file.readline()
	sequence=""

	start_position=(gs_start-1)/50*51+(gs_start-1)%50            # gs_start: 1base

	seq_file.seek(start_position,1)
	seq_length=gs_stop-gs_start+1

	for i in range(0,seq_length):
		letter=seq_file.read(1)
		if letter=="\n":
			letter=seq_file.read(1)
		sequence=sequence+letter

	return sequence

def rev_comp(sequence):
	len_seq=len(sequence)
	rc_seq=""

	for i in range(0,len_seq):
		letter=sequence[len_seq-(i+1)]

		if letter=="A":
			letter="T"
		elif letter=="C":
			letter="G"
		elif letter=="G":
			letter="C"
		elif letter=="T":
			letter="A"
		elif letter=="a":
			letter="t"
		elif letter=="c":
			letter="g"
		elif letter=="g":
			letter="c"
		elif letter=="t":
			letter="a"
		else:
			print("in function rev_comp")
			print("Another letter, exluding ACGTacgt...")
			print sequence
			print letter
			raw_input()
		rc_seq=rc_seq+letter
	return rc_seq


def get_coding_position(gene_exonstart,gene_exonstop,cds_start,cds_stop,breakpoint): # return the SNP position in the gene
	exon_length_1=0
	upper_utr_length_1=0 # genomic position base
	lower_utr_length_1=0 # genomic position base

	exon_breakpoint_pos=0 # genomic position base

	for i in range(0,len(gene_exonstart)):
		exon_length_1+=(gene_exonstop[i]-gene_exonstart[i]+1)
		
		if gene_exonstop[i] < cds_start:
			upper_utr_length_1+= (gene_exonstop[i]-gene_exonstart[i]+1)
		elif gene_exonstart[i] <= cds_start:
			upper_utr_length_1+= ((cds_start-1)-gene_exonstart[i]+1)

		if gene_exonstart[i] > cds_stop:
			lower_utr_length_1 += (gene_exonstop[i]-gene_exonstart[i]+1)
		elif gene_exonstop[i] >= cds_stop:
			lower_utr_length_1 += (gene_exonstop[i] - (cds_stop+1) +1)

		if gene_exonstop[i] < breakpoint:
			exon_breakpoint_pos += (gene_exonstop[i]-gene_exonstart[i]+1)
		elif gene_exonstart[i] <= breakpoint:
			exon_breakpoint_pos += (breakpoint-gene_exonstart[i]+1)

	return [exon_length_1,upper_utr_length_1,lower_utr_length_1,exon_breakpoint_pos]

#def get_prot_seq(sequence,is_donor):
def get_prot_seq(sequence,is_donor,gp_chr):
	cds_length=len(sequence)
	frame=cds_length%3
	frame_num=cds_length/3

	remaining_seq=""
	main_seq=sequence
	protein_seq=""
	
	if cds_length==0:
		return protein_seq

	if is_donor==1:
		if frame!=0:
			remaining_seq=sequence[(-1*frame):].lower()
			main_seq=sequence[0:(-1*frame)]
	else:
		if frame!=0:
			remaining_seq=sequence[0:frame].lower()
			main_seq=sequence[frame:]

	for i in range(0,frame_num):
		this_codon=main_seq[(i*3):(i*3+3)]
		aminoacid=""
		if this_codon in ["TTT","TTC"]:
			aminoacid="F"
		elif this_codon in ["TTA","TTG","CTT","CTC","CTA","CTG"]:
			aminoacid="L"
		elif this_codon in ["ATT","ATC","ATA"]:
			if gp_chr=="MT" and this_codon=="ATA": # Human Mito exception, W. Jia et al., 2008 MBE
				aminoacid="M"
			else:
				aminoacid="I"
		elif this_codon in ["ATG"]:
			aminoacid="M"
		elif this_codon in ["GTT","GTC","GTA","GTG"]:
			aminoacid="V"
		elif this_codon in ["TCT","TCC","TCA","TCG","AGT","AGC"]:
			aminoacid="S"
		elif this_codon in ["CCT","CCC","CCA","CCG"]:
			aminoacid="P"
		elif this_codon in ["ACT","ACC","ACA","ACG"]:
			aminoacid="T"
		elif this_codon in ["GCT","GCC","GCA","GCG"]:
			aminoacid="A"
		elif this_codon in ["TAT","TAC"]:
			aminoacid="Y"
		elif this_codon in ["TAA","TAG","TGA"]:
			if gp_chr=="MT" and this_codon == "TGA":
				aminoacid="W"
			else:
				aminoacid="*"
		elif this_codon in ["CAT","CAC"]:
			aminoacid="H"
		elif this_codon in ["CAA","CAG"]:
			aminoacid="Q"
		elif this_codon in ["AAT","AAC"]:
			aminoacid="N"
		elif this_codon in ["AAA","AAG"]:
			aminoacid="K"
		elif this_codon in ["GAT","GAC"]:
			aminoacid="D"
		elif this_codon in ["GAA","GAG"]:
			aminoacid="E"
		elif this_codon in ["TGT","TGC"]:
			aminoacid="C"
		elif this_codon in ["TGG"]:
			aminoacid="W"
		elif this_codon in ["CGT","CGC","CGA","CGG","AGA","AGG"]:
			if gp_chr=="MT" and this_codon in ["AGA","AGG"]:
				aminoacid="*"
			else:
				aminoacid="R"
		elif this_codon in ["GGT","GGC","GGA","GGG"]:
			aminoacid="G"
		else:
			print("Problem! code 2..no codons...")
			print i
			print "this_codon"
			print this_codon
			print "Main seq:"
			print main_seq
			print "remaining seq:"
			print remaining_seq
			print "all sequence:"
			print sequence
			raw_input()

		protein_seq+=aminoacid
		if aminoacid=="*":
			break

	if is_donor==1:
#		protein_seq=protein_seq+"~"+remaining_seq
		protein_seq=protein_seq
	else:
#		protein_seq=remaining_seq+"~"+protein_seq
		protein_seq=protein_seq

	return protein_seq


#fn="S56_SNP_format.txt"
#fn=raw_input("Input file name: ")

fn=sys.argv[1]

try:
	inputfile=file(fn)
except:
	print("\n\n\t******************** ********************")
	print("\tSNP annotate v0.92 by Young Seok Ju, with MT, last updated Nov 13 2013")
	print("\tjueenome@gmail.com\n")
	print("\tUsage: python 01_SNPannotateXX.py inputfile")
	print("\tTry again...!")
	print("\t******************** ********************")
	sys.exit(0)


version="_annoHg19"
ofn=fn.replace(".txt",version+".txt").replace(".dat",version+".dat").replace(".snp",version+".snp")
if fn==ofn:
	ofn=fn+".annot"
outputfile=file(ofn,"w")

line=inputfile.readline()
while line:
	if "##" in line[0:2]:
		outputfile.write(line)
		line=inputfile.readline()
		continue
	if "#CHROM" in line:
		outputfile.write(line.replace("\n","\tAnnot\tNearest_5p\tNearest_3p\n"))
		line=inputfile.readline() 
	if "\tA\t" in line or "\tC\t" in line or "\tG\t" in line or "\tT\t" in line: #header
		a=1
		break
	else:
		outputfile.write(line.replace("\n","\tAnnot\tNearest_5p\tNearest_3p\n"))
		line=inputfile.readline()
		break 

print("\n\n******************** ********************")
print("SNP annotate v0.91 by Young Seok Ju, with MT, last updated Dec 31 2018")
print("jueenome@gmail.com\n")
print(fn+" is now being annotated using"+ version+"\n")
print("Warning! SNP should be sorted by chromosome then position!")
print("******************** ********************")
print("\n\n")

prev_chr=""
while line:
	print line
	line_split=line.rstrip().split("\t")  # chromosome position ref var (e.g. 1 123213 A G + alpha)
	chr=line_split[0].replace("chr","").replace("23","X").replace("24","Y").replace("MT","M").replace("25","M").replace("M","MT") 
	if ":" in chr:
		chr=chr.split(":")[-1]
	# format: 1, X ...
	pos=int(line_split[1])
	try:
		if prev_chr!=chr:
			prev_chr=chr
			print("Chromosome"+chr)
			infofn="/home/users/ysj/Database/02_genes/01_UCSC_refgene/refGene_chr"+chr+"s.txt"
			infofile=file(infofn)
			infoline=infofile.readline()
	except:
		line=inputfile.readline()
		continue
	wt=line_split[3] # A,C,G,T,N
	var=line_split[4] #A,C,G,T,N
	
	var_length=len(var)-len(wt)
	
	if min(len(var),len(wt))==1: # simple subs or indel
		if var_length>=0: #insertion or subs
			pos1=pos
		else:
			pos1=pos-var_length - 1
	else: #complex
		pos1=pos+len(wt)-1


	gene=""

	anno_info=[] #category,gene,strand,index,context,annotation,blosum90,snperror_flag; 	#category:intergenic,5UTR,3UTR,ncUTR,intronic_noSp,intronic_Sp,synSNP,nsSNP,
											#context: exon #, intron#
											#annotation:for nsSNP V600E(/1200), including total length of protein; others ""
											#flag: if CDS%%3==0 ==> 1, else ==> 4, if wildtype of SNP file != ref seq ==> 2, total_error in protein coding number ==>8 (warning sign! if 4, it may be SNP error!
	snp_info=line.replace("\n","\t")

	flag=0 #0, annotation begins for this line; 1, annotation done for this line

	nearest_gene_5=""
	nearest_gene_3=""
	while infoline:
		infoline_split=infoline.rstrip().split("\t")
#		print(infoline_split)
		infoline_start=int(infoline_split[4])+1
		infoline_stop=int(infoline_split[5])
		gene_index=infoline_split[1]
		gene=infoline_split[12]

		category="intergenic"		

		if flag==0:
			flag=1
			cursor=infofile.tell()
			if cursor > 2000: # 2000, arbitrary large number. byte per line ~ 200, 2000=10 lines
				infofile.seek(-2000,1)
				infoline=infofile.readline()
			else:
				infofile.seek(-1*cursor,1)
			infoline=infofile.readline()
			continue
		if infoline_stop < pos and infoline_stop < pos1:  ######pos1 added for indel
			infoline=infofile.readline()
			nearest_gene_5=gene+"("+gene_index+","+str((pos-infoline_stop+1)/float(1000))+"kb)"			
			continue
		if infoline_start > pos and infoline_start > pos1 :
			nearest_gene_3=gene+"("+gene_index+","+str((pos-infoline_start+1)/float(1000))+"kb)"
			break

		category="intragenic" # intergenic=temporary category

		gene_strand=infoline_split[3]
		if gene_strand!="+" and gene_strand!="-":
			print("****Problem! Strand of gene is not + nor -***")
			print gene_strand
			print infoline
			raw_input()

		gene_cdsstart=int(infoline_split[6])+1
		gene_cdsstop=int(infoline_split[7])
		gene_exoncount=int(infoline_split[8])	
		gene_exonstart=infoline_split[9][0:-1].split(",")
		gene_exonstop=infoline_split[10][0:-1].split(",")
		gene_translationinfo=infoline_split[15][0:-1].split(",") #-1,0,1,2  #-2 for locus information rather than gene

		for i in range(0,gene_exoncount):
			gene_exonstart[i]=int(gene_exonstart[i])+1 #1 base
			gene_exonstop[i]=int(gene_exonstop[i])
			gene_translationinfo[i]=int(gene_translationinfo[i])

		#***category testing!
		category_flag="intronic_noSp"
		for i in range(0,gene_exoncount):
			if pos>=gene_exonstart[i] and pos<=gene_exonstop[i]:
				if gene_translationinfo.count(-1)==gene_exoncount: # noncoding gene
					category_flag="ncUTR"
				elif -2 in gene_translationinfo:
					category_flag="LocInfo" # locus info... Such as attachment site... etc
				elif pos < gene_cdsstart:
					if gene_strand=="+":
						category_flag="5UTR"
					else:
						category_flag="3UTR"
				elif pos > gene_cdsstop:
					if gene_strand=="+":
						category_flag="3UTR"
					else:
						category_flag="5UTR"
				else:
						category_flag="CDS" # CDS; temporary category
			elif (pos >= (gene_exonstart[i]-2) and pos <= (gene_exonstart[i]-1)) or (pos >= (gene_exonstop[i]+1) and pos <= (gene_exonstop[i]+2)):
						category_flag="intronic_Sp"

		category_flag1="intronic_noSp" #####################added for indel
		for i in range(0,gene_exoncount):
			if pos1>=gene_exonstart[i] and pos1<=gene_exonstop[i]:
				if gene_translationinfo.count(-1)==gene_exoncount: # noncoding gene
					category_flag1="ncUTR"
				elif -2 in gene_translationinfo:
					category_flag1="LocInfo" # locus info... Such as attachment site... etc
				elif pos1 < gene_cdsstart:
					if gene_strand=="+":
						category_flag1="5UTR"
					else:
						category_flag1="3UTR"
				elif pos1 > gene_cdsstop:
					if gene_strand=="+":
						category_flag1="3UTR"
					else:
						category_flag1="5UTR"
				else:
						category_flag1="CDS" # CDS; temporary category
			elif (pos1 >= (gene_exonstart[i]-2) and pos1 <= (gene_exonstart[i]-1)) or (pos1 >= (gene_exonstop[i]+1) and pos1 <= (gene_exonstop[i]+2)):
						category_flag1="intronic_Sp"

		category_membership=0
		pos2=0
		category_flag_final=""
		if category_flag1==category_flag: #pos, pos1 in the same category
			category_membership=0
			category_flag_final=category_flag
			pos2=pos
		else:
			category_membership=1
			if category_flag=="intronic_Sp":
				pos2=pos
				category_flag_final=category_flag
			elif category_flag1=="intronic_Sp":
				pos2=pos1
				category_flag_final=category_flag1
			elif category_flag=="CDS":
				pos2=pos
				category_flag_final="intronic_Sp"
			elif category_flag1=="CDS":
				pos2=pos1
				category_flag_final="intronic_Sp"
			else:
				pos2=pos # any will do
				category_flag_final=category_flag

		#coding_positions=get_coding_position(gene_exonstart,gene_exonstop,gene_cdsstart,gene_cdsstop,pos)
		coding_positions=get_coding_position(gene_exonstart,gene_exonstop,gene_cdsstart,gene_cdsstop,pos2)

		exon_length=coding_positions[0]
		up_utrlength=coding_positions[1]
		down_utrlength=coding_positions[2]
		cds_length=exon_length-up_utrlength-down_utrlength
		breakpoint_length=coding_positions[3]

		snperror_flag=1
		aminoacid_length=cds_length/3

		if cds_length%3 !=0:
			print "***Warning! CDS length is not 3n!***"
			print(infoline.rstrip())
			print("Exon length:"+str(exon_length)+"bp")
#			print up_utrlength
#			print down_utrlength
			print("CDS length:"+str(cds_length)+"bp")
			snperror_flag+=4

		#******Exon, intron number ****************#
		print pos
		print pos1
		print pos2
		print category_flag
		print category_flag1
		print category_flag_final
		exon_num=gene_exoncount-1
		intron_flag=0 # 1 intron
		for i in range(0,gene_exoncount):
			if (int(gene_exonstart[i])+1-pos2)*(int(gene_exonstop[i])-pos2) <=0:
				exon_num=i
				break
			if pos2 < (int(gene_exonstart[i])+1):
				exon_num=i
				intron_flag = 1
				break

		if intron_flag==1: #intron
			if gene_strand=="+":
				position_number="intron"+str(exon_num)
			elif gene_strand=="-":
				position_number="intron"+str(gene_exoncount-exon_num)
			else:
				print "problem! Intron number...!"
				print gene_strand
				print line
				print infoline
				raw_input()
			position_start=int(gene_exonstop[exon_num-1])+1
			position_stop=int(gene_exonstart[exon_num])+1-1
		else:
			if gene_strand=="+":
				if category_flag_final=="LocInfo":
					position_number="segment"+str(exon_num+1)
				else:
					position_number="exon"+str(exon_num+1)
			elif gene_strand=="-":
				if category_flag_final=="LocInfo":
					position_number="segment"+str(exon_num+1)
				else:
					position_number="exon"+str(gene_exoncount-exon_num)
			else:
				print "problem! Exon number"
				print gene_strand
				print line
				print infoline
				raw_input()
			print exon_num
			print infoline_split
			position_start=int(gene_exonstart[exon_num])+1
			position_stop=int(gene_exonstop[exon_num])
	
#			outputfile.write(gene+"("+gene_index+","+gene_strand+"strand),"+position_number+"("+chr+":"+str(position_start)+"-"+str(position_stop)+";")

		if ("UTR" in category_flag_final) or ("intronic" in category_flag_final) or ("LocInfo" in category_flag_final):
#			anno_info=[] #category,gene,index,context,annotation,blosum,snperror_flag; 	#category:intergenic,5UTR,3UTR,ncUTR,intronic_noSp,intronic_Sp,synSNP,nsSNP,
			category=category_flag_final
			if gene_index in ["MT-TF","MT-TV","MT-TL1","MT-TI","MT-TQ","MT-TM","MT-TW","MT-TA","MT-TN","MT-TC","MT-TY","MT-TS1","MT-TD","MT-TK","MT-TG","MT-TR","MT-TH","MT-TS2","MT-TL2","MT-TE","MT-TT","MT-TP"]: # for MT-tRNA annotation
				tRNAinfo=tRNAlocus(gene_index,gene_strand,pos) ####FIXING
				anno_info.append(category+","+gene+","+gene_strand+","+gene_index+","+position_number+","+tRNAinfo[0]+","+tRNAinfo[1]+","+str(snperror_flag)+";")
			else:
				anno_info.append(category+","+gene+","+gene_strand+","+gene_index+","+position_number+",,,"+str(snperror_flag)+";")
			infoline=infofile.readline()
			continue

		# CDS snp (synonymous or non-synonymous

		#***************get sequence************************
		raw_sequence=""
		for i in range(0, gene_exoncount):
			gs_start1=gene_exonstart[i]
			gs_stop1=gene_exonstop[i]
			raw_sequence=raw_sequence+get_sequence(chr,gs_start1,gs_stop1)

		RNAseq=""
		CDSseq=""
		PROTEINseq=""
	
		WARNINGmessage=""

		CDSseq=raw_sequence[up_utrlength:(down_utrlength*-1)].upper()
#		CDSseq_var=raw_sequence[up_utrlength:(breakpoint_length-1)].upper()+var.upper()+raw_sequence[(breakpoint_length):(down_utrlength*-1)].upper()
		CDSseq_var=raw_sequence[up_utrlength:(breakpoint_length-1)].upper()+var.upper()+raw_sequence[(breakpoint_length+len(wt)-1):(down_utrlength*-1)].upper()

		RNAseq=raw_sequence[0:up_utrlength].lower()+CDSseq+raw_sequence[(down_utrlength*-1):].lower()
		RNAseq_var=raw_sequence[0:up_utrlength].lower()+CDSseq_var+raw_sequence[(down_utrlength*-1):].lower()
		if down_utrlength==0:
			CDSseq=raw_sequence[up_utrlength:].upper()
			CDSseq_var=raw_sequence[up_utrlength:(breakpoint_length-1)].upper()+var.upper()+raw_sequence[(breakpoint_length+len(wt)-1):].upper()
			RNAseq=raw_sequence[0:up_utrlength].lower()+CDSseq
			RNAseq_var=raw_sequence[0:up_utrlength].lower()+CDSseq_var

		print raw_sequence
		if raw_sequence[breakpoint_length-1].upper() !=	wt:
			print("****Warning! wildtype in human reference genome and wt allele in your file is different!*****")
			#:print line
			#print gene
			#print wt
			#print("raw_sequence[breakpoint_length-1]")
			#print raw_sequence[breakpoint_length-1]
			#print breakpoint_length
			#print("raw_sequence")
			#print raw_sequence
			#print raw_sequence[741]
			#print raw_sequence[952]
			#print raw_sequence[953]
			#print raw_sequence[954]
			#print raw_sequence[955]
			#print raw_sequence[956]
			#print len(raw_sequence)
			snperror_flag+=2		

#		print CDSseq
#		print CDSseq_var
#		print "111"
#		k=raw_input()



		if gene_strand=="-":
			CDSseq=rev_comp(CDSseq)
			RNAseq=rev_comp(RNAseq)
			CDSseq_var=rev_comp(CDSseq_var)
			RNAseq=rev_comp(RNAseq_var)

		PROTEINseq=get_prot_seq(CDSseq,1,chr)
		PROTEINseq_var=get_prot_seq(CDSseq_var,1,chr)

		if PROTEINseq[-1]!="*":
			PROTEINseq=PROTEINseq+"*"
		if PROTEINseq_var[-1]!="*":
			PROTEINseq_var=PROTEINseq+"*"

#		PROTEINlength=len(CDSseq)/3-1 #exclude stop codon
		PROTEINlength=len(PROTEINseq)-1

		if len(wt)==1 and len(var)==1: #substitutions		
			AA_altered=(breakpoint_length-up_utrlength+2)/3
			if gene_strand=="-":
				AA_altered=PROTEINlength-AA_altered+2

			try:
				AA_wt=PROTEINseq[AA_altered-1]
				AA_var=PROTEINseq_var[AA_altered-1]
			except:
				AA_wt="*"
				AA_var="*"
				snperror_flag+=8

			if AA_wt==AA_var:
				category_flag="synSNP"
			else:
				category_flag="nsSNP"
			annotation=AA_wt+str(AA_altered)+AA_var
			annotation1=PROTEINseq[max((AA_altered-10),0):(AA_altered-1)].lower()+"["+AA_wt+">"+AA_var+"]"+PROTEINseq[AA_altered:min((PROTEINlength),AA_altered+9)].lower()
			blosum_score=blosum90(AA_wt,AA_var)	
		else: #complex or indels
			if PROTEINseq==PROTEINseq_var:
				category_flag="syn_multiplesubs"
				annotation="NA"
				annotation1="NA"
			else:
				if (var_length)%3==0:
					category_flag="inframe_CDS_INDEL"
				else:
					category_flag="frameshift_CDS_INDEL"
				identical_aa_num=0
				identical_aa_num_backward = 0
				for amino_acid_num in range(0,PROTEINlength):
					if PROTEINseq[amino_acid_num]==PROTEINseq_var[amino_acid_num]:
						identical_aa_num+=1
						continue
					break
				for amino_acid_num in range(1,PROTEINlength+2):
					if PROTEINseq[amino_acid_num*-1]==PROTEINseq_var[amino_acid_num*-1]:
						identical_aa_num_backward -= 1 
						continue
					break
				print category_flag
				print infoline_split
				print PROTEINseq
				print PROTEINseq_var
				if identical_aa_num_backward!=0:
					annotation=str(identical_aa_num+1)+PROTEINseq[identical_aa_num:identical_aa_num_backward]+">"+PROTEINseq_var[identical_aa_num:identical_aa_num_backward]

					#HOMOPOLYMER "ASV>ASVASV"
					if PROTEINseq[identical_aa_num:identical_aa_num_backward]=="" and PROTEINseq_var[identical_aa_num:identical_aa_num_backward]=="":
					#	print "UES"
						identical_aa_num_backward=identical_aa_num_backward+len(PROTEINseq_var)-len(PROTEINseq)
						annotation=str(identical_aa_num+1)+PROTEINseq[identical_aa_num:identical_aa_num_backward]+">"+PROTEINseq_var[identical_aa_num:identical_aa_num_backward]
				else:
					annotation=str(identical_aa_num+1)+PROTEINseq[identical_aa_num:]+">"+PROTEINseq_var[identical_aa_num:]

			
				if identical_aa_num_backward<-9:
					annotation1=PROTEINseq[max(identical_aa_num-9,0):identical_aa_num].lower()+"["+annotation.upper()+"]"+PROTEINseq[identical_aa_num_backward:(identical_aa_num_backward+9)].lower()
				else :# identical_aa_num_backward<-1:
					annotation1=PROTEINseq[max(identical_aa_num-9,0):identical_aa_num].lower()+"["+annotation.upper()+"]"+PROTEINseq[identical_aa_num_backward:].lower()
#				else:a
#					annotation1=PROTEINseq[max(identical_aa_num-9,0):identical_aa_num].lower()+annotation.upper()
#				print PROTEINseq
#				print PROTEINseq_var
#				print identical_aa_num
#				print identical_aa_num_backward
#				print annotation
	#			print annotation1
#				raw_input("!!!indel protein sequences...")
			blosum_score=-99
			
#		anno_info=[] #category,gene,index,context,annotation,blosum,snperror_flag; 	#category:intergenic,5UTR,3UTR,ncUTR,intronic_noSp,intronic_Sp,synSNP,nsSNP,
		category=category_flag
		anno_info.append(category+","+gene+","+gene_strand+","+gene_index+","+position_number+","+annotation+","+annotation1+","+str(blosum_score)+","+str(snperror_flag)+";")

		infoline=infofile.readline()

	if anno_info==[]:
		anno_info.append("intergenic,,,,,,,;")

	outputfile.write(snp_info+"".join(anno_info)+"\t"+nearest_gene_5+"\t"+nearest_gene_3+"\n")
	line=inputfile.readline()



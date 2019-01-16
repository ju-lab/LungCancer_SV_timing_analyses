#This script was made 2018-06-29: modified from snv_timing_probability.py which use the definition of earlySNV as mutCN==majCN
#In this script, the definition of earlySNV is mutCN >= round(majCN*3/4,0)
#181231 input column edit

import sys,math
import scipy.stats as ss
print('###'+sys.argv[0])
def calc_binomialP(dp, var, tCN, MCN, purity):
	#dp, var, tCN, MCN should be integer
	dp=int(dp);var=int(var);tCN=int(tCN);MCN=int(MCN)
	purity=float(purity)
	tfraction=tCN*purity/float(tCN*purity + 2*(1-purity))
	nE=0;nL=0;nS=0
	if tCN==0:
		return ['.','.','.','.']

	cutoff=int(math.ceil(MCN*3/float(4)))
	for n in range(cutoff, MCN+1):
		hh=ss.binom(dp, tfraction*n/tCN)
		if dp*tfraction*n/tCN < var:
			for r in range(var, dp+1):
				nE += hh.pmf(r)
		elif dp*tfraction*n/tCN >= var:
			for r in range(0, var+1):
				nE += hh.pmf(r)

	for n in range(1, cutoff):
		hh=ss.binom(dp, tfraction*n/tCN)
		if dp*tfraction*n/tCN < var:
			for r in range(var, dp+1):
				nL += hh.pmf(r)
		elif dp*tfraction*n/tCN >= var:
			for r in range(0, var+1):
				nL += hh.pmf(r)
	hh=ss.binom(dp,tfraction*0.5/tCN)
	if dp*tfraction*0.5/tCN < var:
		for r in range(var, dp+1):
			nS += hh.pmf(r)
	elif dp*tfraction*0.5/tCN >= var:
		for r in range(0, var+1):
			nS += hh.pmf(r)
	totaln=nE+nL+nS
	pC=(nE+nL)/float(totaln);pE=nE/float(totaln); pL=nL/float(totaln); pS=nS/float(totaln)
	if nL == 0:
		return [str(pC),'.','.',str(pS)]
	else:
		return [str(pC), str(pE), str(pL), str(pS)]

	

print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.proba','w')
purity=float(sys.argv[2])
in_line=in_file.readline().strip()
while in_line:
	if in_line[0:6]=='#CHROM':
		out_file.write(in_line+'\tpClonal\tpEarly\tpLate/minor\tpSubclonal\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		refn=in_indi[8]
		varn=in_indi[9]
		tCN=in_indi[10]
		MCN=in_indi[11]
		if MCN=='NA':
			out_file.write(in_line+'\t.\t.\t.\t.\n')
			in_line=in_file.readline().strip()
			continue
		dp= int(refn)+int(varn)
		var=int(varn)
		p_list=calc_binomialP(dp, var, tCN, MCN, purity)
		out_file.write(in_line+'\t'+'\t'.join(p_list)+'\n')
	in_line=in_file.readline().strip()


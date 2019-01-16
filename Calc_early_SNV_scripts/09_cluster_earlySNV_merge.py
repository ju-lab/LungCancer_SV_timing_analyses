#Arg1: cluster cl_time

import sys, math
print('###'+sys.argv[0])
ampsizeco=10000000
clsvn_co=5
print(sys.argv[1])
in_file=open(sys.argv[1])  # cl_time
out_file=open(sys.argv[1]+'.merge','w')
header_list=['#cluster', 'chr_list', 'dist', 'early_SNV', 'early_SNV_range', 'early_SNVrate/Mb', 'early_SNVrate_range', 'early_SNV_var','BPnum', 'cluster_SVnum']
out_file.write('\t'.join(header_list)+'\n')
cl_list=[]
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':'blank'
	else:
		in_indi=in_line.split('\t')
		cl=in_indi[0]
		cl_list.append(cl)
	in_line=in_file.readline().strip()

cl_list=list(set(cl_list))

for clst in cl_list:
	print(clst)
	target_chr_list=[]
	dist_sum=0; earlyn_sum=0;earlyv_sum=0;bp_sum=0
	in_file.seek(0)
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':'blank'
		else:
			in_indi=in_line.split('\t')
			cl=in_indi[0];chr1=in_indi[1];majCN=int(in_indi[5]);minCN=int(in_indi[6])
			dist=int(in_indi[8]);broadinfo=in_indi[7]
			if cl==clst and majCN >=2 and broadinfo=='.' :
				target_chr_list.append(chr1)
				earlyn=float(in_indi[9]); earlyvar=float(in_indi[13])
				bpnum=int(in_indi[14]); cl_svn=int(in_indi[15])
				if minCN <= majCN*0.75:
					dist_sum += dist; earlyn_sum +=2*earlyn; earlyv_sum += 2*earlyvar; bp_sum += bpnum
				else:
					dist_sum += dist; earlyn_sum +=earlyn; earlyv_sum += earlyvar ; bp_sum += bpnum
					
		in_line=in_file.readline().strip()
	target_chr_list=list(set(target_chr_list))
	print(dist_sum)
	if dist_sum >= ampsizeco and cl_svn >=cl_svn:
		earlyr=round(earlyn_sum*1000000/float(dist_sum),4)
		learlyn_sum = earlyn_sum - 1.96*math.sqrt(earlyv_sum)
		hearlyn_sum = earlyn_sum + 1.96*math.sqrt(earlyv_sum)
		learlyr=round(learlyn_sum*1000000/float(dist_sum),4)
		hearlyr=round(hearlyn_sum*1000000/float(dist_sum),4)
		earlyar=round(earlyr/float(tclr),4)
		learlyar=round(learlyr/float(tclr),4)
		hearlyar=round(hearlyr/float(tclr),4)
		earlyar2=round(earlyr/float(tsnvr),4)
		learlyar2=round(learlyr/float(tsnvr),4)
		hearlyar2=round(hearlyr/float(tsnvr),4)
		chr_list=';'.join(target_chr_list)

		info_list=[clst,chr_list,str(dist_sum), str(earlyn_sum), str(learlyn_sum)+':'+str(hearlyn_sum),str(earlyr), str(learlyr)+':'+str(hearlyr),str(tcln),str(tclr),str(tsnvn),str(tsnvr),str(earlyar),str(learlyar)+':'+str(hearlyar),str(earlyar2),str(learlyar2)+':'+str(hearlyar2),tclnr, tclnrr, str(earlyv_sum), str(bp_sum),str(cl_svn)]
		out_file.write('\t'.join(info_list)+'\n')

#Arg1: earlySNV  (ampseg.earlysnv.txt)
#Arg2: sampleid
#Arg3: outDit

import sys, math
from operator import itemgetter
print('###'+sys.argv[0])
print(sys.argv[1])
seg_file=open(sys.argv[1]) # segment with time info, *.ampseg.timeci.txt
out_file=open(sys.argv[3]+'/'+sys.argv[2]+'.broad_earlysnv','w') # id
header_list=['#broad_amp','dist','early_SNV','early_SNV_range', 'early_SNVrate/Mb', 'early_SNVrate_range', 'early_SNV_var']
out_file.write('\t'.join(header_list)+'\n')
seg_file.seek(0)
seg_line=seg_file.readline().strip()
broad_list=[]
while seg_line:
	if seg_line[0]=='#':
		'blank'
	else:
		seg_indi=seg_line.split('\t')
		broadinfo=seg_indi[6]
		if broadinfo != '.':
			broad_chr=int((((broadinfo.split('p')[0]).split('q')[0]).split('W')[0].replace('X','23')).replace('Y','24'))
			if 'p' in broadinfo:
				broad_idx=0
			elif 'q' in broadinfo:
				broad_idx=1
			else:
				broad_idx=2
			if [broad_chr,broad_idx,broadinfo] not in broad_list:
				broad_list.append([broad_chr,broad_idx,broadinfo])
	seg_line=seg_file.readline().strip()

broad_list.sort(key=itemgetter(0,1))

for [chr1,idx,amp] in broad_list:
	dist_sum=0; earlyn_sum=0;earlyV_sum=0
	seg_file.seek(0)
	seg_line=seg_file.readline().strip()
	while seg_line:
		if seg_line[0]=='#':
			'blank'
		else:
			seg_indi=seg_line.split('\t')
			broadinfo=seg_indi[6]; majCN=int(seg_indi[4]);minCN=int(seg_indi[5])
			if broadinfo==amp and majCN>=2:
				dist=int(seg_indi[7]);earlyn=float(seg_indi[8]); earlyV=float(seg_indi[12])
				dist_sum += dist; 
				if minCN <= majCN*0.75:
					earlyn_sum += 2*earlyn; earlyV_sum += 2*earlyV
				else:
					earlyn_sum += earlyn; earlyV_sum += earlyV
		seg_line=seg_file.readline().strip()

	earlyr=round(earlyn_sum*1000000/float(dist_sum),4)
	learlyn_sum = earlyn_sum - 1.96*math.sqrt(earlyV_sum)
	hearlyn_sum = earlyn_sum + 1.96*math.sqrt(earlyV_sum)
	learlyr=round(learlyn_sum*1000000/float(dist_sum),4)
	hearlyr=round(hearlyn_sum*1000000/float(dist_sum),4)
	info_list=[amp,str(dist_sum), str(earlyn_sum), str(learlyn_sum)+':'+str(hearlyn_sum),str(earlyr), str(learlyr)+':'+str(hearlyr),str(earlyV_sum)]
	if 'X' not in amp and 'Y' not in amp:
		out_file.write('\t'.join(info_list)+'\n')

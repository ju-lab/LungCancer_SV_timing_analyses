#Args1  5Mclustered file

import sys
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.edit','w')
header_list=['#CHR1', 'POS1', 'CHR2', 'POS2', 'mh', 'terinfo', 'svtype', 'precise', 'class1', 'conc1', 'allfrag', 'spfrag', 'safrag', 'ref1', 'ref2', 'mapq1', 'mapq2', 'tbginfo1', 'tbginfo2', 'nbginfo2', 'nbginfo2', 'pos_dist1', 'pos_dist2', 'dellyid', 'unknown_bp1', 'unknown_bp2', 'cluster']
out_file.write('\t'.join(header_list)+'\n')

in_line=in_file.readline().strip()
while in_line:
	if in_line[0:5]=='#clus':
		cluster_num=in_line.split('cluster')[1]
	else:
		out_file.write(in_line+'\t'+cluster_num+'\n')
	in_line=in_file.readline().strip()
	

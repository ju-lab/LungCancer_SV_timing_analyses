import sys

sv_dist=1000
near_dist=10000
bal_dist=500
overlap=1

overlap=(-1)*overlap+1

in_file=open(sys.argv[1])
print(sys.argv[1])
out_file=open(sys.argv[1]+'.split_simple','w')

in_dic={}

in_line=in_file.readline().strip()
c_list=[]
n=0;m=0
while in_line:
	if in_line[0:6]=='#clust':
		m+=1
		cluster_num=in_line.split('cluster')[1]
		c_list.append(int(cluster_num))
		in_dic[cluster_num]=[]
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		n+=1
		in_dic[cluster_num].append(in_line)
	in_line=in_file.readline().strip()
print(n)
print(m)

new_cluster_list=[]
for clust in in_dic.keys():
	for line1 in in_dic[clust]:
		chr1=line1.split('\t')[0]
		pos1=line1.split('\t')[1]
		chr2=line1.split('\t')[2]
		pos2=line1.split('\t')[3]
		ter1=line1.split('\t')[5][0]
		ter2=line1.split('\t')[5][-1]
		svtype1=line1.split('\t')[6]
		idx1='\t'.join([chr1,pos1,chr2,pos2,ter1,ter2])
		adj_sv1=0;adj_sv2=0
		simple_del='F';simple_dup='F';simple_inv='F';simple_tra='F';simple_event='F'
		#check nearby SV
		for line2 in in_dic[clust]:
			chr3=line2.split('\t')[0]
			pos3=line2.split('\t')[1]
			chr4=line2.split('\t')[2]
			pos4=line2.split('\t')[3]
			ter3=line2.split('\t')[5][0]
			ter4=line2.split('\t')[5][-1]
			idx2='\t'.join([chr3,pos3,chr4,pos4,ter3,ter4])
			if idx1 == idx2:
				continue
			if chr1==chr3 and abs(int(pos1)-int(pos3)) < near_dist:
				adj_sv1+=1
			if chr4 !='.':
				if chr1==chr4 and abs(int(pos1)-int(pos4)) < near_dist:
					adj_sv1+=1
			if chr2 != '.':
				if chr2 == chr3 and abs(int(pos2)-int(pos3)) < near_dist:
					adj_sv2+=1
				if chr2 == chr4 and abs(int(pos2) -int(pos4)) < near_dist:
					adj_sv2+=1
		if adj_sv1 >= 2 or adj_sv2 >=2:
			continue
		#check simple_sv
		if svtype1 == 'DEL' and int(pos2)-int(pos1) < sv_dist and adj_sv1==0 and adj_sv2 ==0:
			new_cluster_list.append([])
			new_cluster_list[-1].append(line1)
			del in_dic[clust][in_dic[clust].index(line1)]
		elif svtype1=='DUP' and int(pos2)-int(pos1) < sv_dist and adj_sv1==0 and adj_sv2==0:
			new_cluster_list.append([])
			new_cluster_list[-1].append(line1)
			del in_dic[clust][in_dic[clust].index(line1)]
		elif svtype1=='INV' and ter1=='3' and ter2=='3':
			for line2 in in_dic[clust]:
				chr3=line2.split('\t')[0]
				pos3=line2.split('\t')[1]
				chr4=line2.split('\t')[2]
				pos4=line2.split('\t')[3]
				ter3=line2.split('\t')[5][0]
				ter4=line2.split('\t')[5][-1]
				svtype2=line2.split('\t')[6]
				idx2='\t'.join([chr3,pos3,chr4,pos4,ter3,ter4])
				if idx1 == idx2:
					continue
				if svtype2 == 'INV' and ter3=='5' and ter4 =='5' and chr1==chr3 and chr2==chr4 and int(pos1)-int(pos3) >=overlap and int(pos1)-int(pos3) < bal_dist and int(pos2)-int(pos4) >=overlap and int(pos2)-int(pos4) < bal_dist and adj_sv1 == 1 and adj_sv2==1:
					new_cluster_list.append([])
					new_cluster_list[-1].append(line1)
					new_cluster_list[-1].append(line2)
					del in_dic[clust][in_dic[clust].index(line1)]
					del in_dic[clust][in_dic[clust].index(line2)]
		elif svtype1=='INV' and ter1=='5' and ter2=='5':
			for line2 in in_dic[clust]:
				chr3=line2.split('\t')[0]
				pos3=line2.split('\t')[1]
				chr4=line2.split('\t')[2]
				pos4=line2.split('\t')[3]
				ter3=line2.split('\t')[5][0]
				ter4=line2.split('\t')[5][-1]
				svtype2=line2.split('\t')[6]
				idx2='\t'.join([chr3,pos3,chr4,pos4,ter3,ter4])
				if idx1 == idx2:
					continue
				if svtype2 == 'INV' and ter3=='3' and ter4=='3' and chr1==chr3 and chr2==chr4 and int(pos3)-int(pos1) >=overlap and int(pos3)-int(pos1) < bal_dist and int(pos4)-int(pos2) >=overlap and int(pos4)-int(pos2) < bal_dist and adj_sv1==1 and adj_sv2==1:
					new_cluster_list.append([])
					new_cluster_list[-1].append(line1)
					new_cluster_list[-1].append(line2)
					del in_dic[clust][in_dic[clust].index(line1)]
					del in_dic[clust][in_dic[clust].index(line2)]
		elif svtype1=='TRA':
			if ter1=='3':
				pos1n=int(pos1)*(-1)
			elif ter1=='5':
				pos1n=int(pos1)
			if ter2=='3':
				pos2n=int(pos2)*(-1)
			elif ter2=='5':
				pos2n=int(pos2)
			for line2 in in_dic[clust]:
				chr3=line2.split('\t')[0]
				pos3=line2.split('\t')[1]
				chr4=line2.split('\t')[2]
				pos4=line2.split('\t')[3]
				ter3=line2.split('\t')[5][0]
				ter4=line2.split('\t')[5][-1]
				svtype2=line2.split('\t')[6]
				idx2='\t'.join([chr3,pos3,chr4,pos4,ter3,ter4])
				if idx1 == idx2:
					continue
				if svtype2 == 'TRA':
					if ter3=='3':
						pos3n=int(pos3)*(-1)
					elif ter3=='5':
						pos3n=int(pos3)
					if ter4=='3':
						pos4n=int(pos4)*(-1)
					elif ter4=='5':
						pos4n=int(pos4)
					if chr1==chr3 and chr2==chr4 and pos1n+pos3n >= overlap and pos1n+pos3n < bal_dist and pos2n+pos4n >=0 and pos2n+pos4n < bal_dist:
							new_cluster_list.append([])
							new_cluster_list[-1].append(line1)
							new_cluster_list[-1].append(line2)
							del in_dic[clust][in_dic[clust].index(line1)]
							del in_dic[clust][in_dic[clust].index(line2)]
					elif chr1==chr4 and chr2==chr3 and pos1n+pos4n >= overlap and pos1n+pos4n < bal_dist and pos2n+pos3n >=0 and pos2n+pos3n < bal_dist:
							new_cluster_list.append([])
							new_cluster_list[-1].append(line1)
							new_cluster_list[-1].append(line2)
							del in_dic[clust][in_dic[clust].index(line1)]
							del in_dic[clust][in_dic[clust].index(line2)]

n=0;m=0
for clust in in_dic.keys():
	if len(in_dic[clust]) > 0:
		m+=1
		out_file.write('#cluster'+str(m)+'\n')
		for line in in_dic[clust]:
			n+=1
			out_file.write(line+'\n')

for clust in new_cluster_list:
	m+=1
	out_file.write('#cluster'+str(m)+'\n')
	for line in clust:
		n+=1
		out_file.write(line+'\n')

print(n)
print(m)

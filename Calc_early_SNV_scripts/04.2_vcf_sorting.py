import sys,os
from operator import itemgetter
print('###'+sys.argv[0])
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.sorted','w')
in_line=in_file.readline().strip()
input_list=[]
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		if in_indi[0][0:2]=='GL':
			chr1=float(in_indi[0][2:])+25
		else:
			chr1=int(((in_indi[0].replace('X','23')).replace('Y','24')).replace('MT','25'))
		input_list.append([chr1,int(in_indi[1]),in_line])
	in_line=in_file.readline().strip()


input_list.sort(key=itemgetter(0,1))

prev_chr=0
prev_pos=0
prev_line='blank'
for [a,b,c] in input_list:
	if prev_chr==a and prev_pos==b:
		print('#Warning repeated position# ')
		print(prev_line)
		print(c)
	out_file.write(c+'\n')
	prev_chr=a
	prev_pos=b
	prev_line=c


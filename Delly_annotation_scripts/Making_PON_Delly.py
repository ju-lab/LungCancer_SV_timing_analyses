import sys,argparse
from operator import itemgetter
parser = argparse.ArgumentParser()
parser.add_argument("input", help= "ID1\\t/path/to/dellyoutput1.vcf\\nID2\\t/path/to/dellyoutput2.vcf\\n....")
parser.add_argument("-n", "--name", default="PON.delly.txt", help='Name of ouput file')
parser.add_argument("-c", "--column", type=int, default='11', help="Number of column with Normal information in input Delly files")
args=parser.parse_args()

fn_file=open(args.input)
fn_line=fn_file.readline().strip()
out_file=open(args.name,'w')
ncol=int(args.column)-1
out_file.write('#CHROM\tPOS\tID\tREF\tALT\tFILTER\tQUAL\tINFO\tFORMAT\tNormal\tSampleID\n')
dt=[]
chr_pos_dic={}
n=0
while fn_line:
	f_id=fn_line.split('\t')[0]
	f_path=fn_line.split('\t')[1]
	in_file=open(f_path)
	in_line=in_file.readline().strip()
	while in_line:
		in_indi=in_line.split('\t')
		if in_line[0:6]=='#CHROM':
			print(in_indi[ncol])
		elif in_line[0]=='#':
			'blank'
		else:
			n_dv=int(in_indi[ncol].split(':')[9])
			n_rv=int(in_indi[ncol].split(':')[11])
			n_chr2=(in_indi[7].split('CHR2=')[1]).split(';')[0]
			n_end=(in_indi[7].split('END=')[1]).split(';')[0]
			ter1=(in_indi[7].split('CT=')[1]).split('to')[0]
			ter2=((in_indi[7].split('CT=')[1]).split('to')[1]).split(';')[0]
			svtype=in_indi[2][0:3]
				
			if in_indi[0][0:2]=='GL':
				chr=25+float(in_indi[0][2:])
			elif in_indi[0][0:3]=='NC_':
				chr=26+(float(in_indi[0][3:])/1000000)
			elif in_indi[0]=='hs37d5':
				chr = 27
			else:
				chr=int(((in_indi[0].replace('X','23')).replace('Y','24')).replace('MT','25'))

			if n_chr2[0:2]=='GL':
				chr2=25+float(n_chr2[2:])
			elif n_chr2[0:3]=='NC_':
				chr2=26+(float(n_chr2[3:])/1000000)
			elif n_chr2 == 'hs37d5':
				chr2 = 27
			else:
				chr2=int(((n_chr2.replace('X','23')).replace('Y','24')).replace('MT','25'))
			org_info='\t'.join(in_indi[0:9])+'\t'+in_indi[ncol]+'\t'+f_id
			rvs_info=n_chr2+'\t'+n_end+'\t'+in_indi[2]+'\t.\t'+'\t'.join(in_indi[4:7])+'\t'+in_indi[7].split('CHR2=')[0]+'CHR2='+in_indi[0]+';END='+in_indi[1]+';PE='+(in_indi[7].split(';PE=')[1]).split(';CT=')[0]+';CT='+ter2+'to'+ter1+';'+';'.join((in_indi[7].split(';CT=')[1]).split(';')[1:])+'\t'+in_indi[8]+'\t'+in_indi[ncol]+'\t'+f_id
			if n_dv + n_rv > 0:    # Cut off for enlisting
				if svtype != 'TRA':
					dt.append([chr,int(in_indi[1]),chr2, int(n_end),org_info])   # Append original info
				elif svtype == 'TRA':
					if chr < chr2:
						dt.append([chr,int(in_indi[1]),chr2, int(n_end),org_info])   # Append original info
					elif chr > chr2:
						dt.append([chr2,int(n_end),chr, int(in_indi[1]),rvs_info]) # Append reverse info
		in_line=in_file.readline().strip()
	fn_line=fn_file.readline().strip()
				
print('sorting')
dt.sort(key=itemgetter(0,1,2,3))

print('writing')
for [a,b,c,d,e] in dt:
	out_file.write(e+'\n')



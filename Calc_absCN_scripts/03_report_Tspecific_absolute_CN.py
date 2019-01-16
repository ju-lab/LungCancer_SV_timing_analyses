#Arg1: tumor 100kbcov
#Arg2: normal 100kbcov
#Arg3: tumor coverage
#Arg4: normal coverage
#Arg5: aberrant cellular fraction
#Arg6: ploidy

import sys

def absolute_copy(cov1,cov2):
		if cov2 < 0.1 :
				return -1
		CN=(cov1/C_cov/cov2*N_cov-(1-acf))*ploidy/acf
		return CN

print("Tumorfile\tNormalfile\tTumorCov\tNormalCov\tACF\tploidy\n")

fnC=sys.argv[1]
fnN=sys.argv[2]
C_cov=float(sys.argv[3]) #43.33
N_cov=float(sys.argv[4]) #22.70
acf=float(sys.argv[5])
ploidy=float(sys.argv[6])

inputfile1=file(fnC)
line1=inputfile1.readline() #header
line1=inputfile1.readline() #header
line1=inputfile1.readline() 

inputfile2=file(fnN)
line2=inputfile2.readline() #header
line2=inputfile2.readline() #header
line2=inputfile2.readline() 

ofn=fnC+".absCN"
outputfile=file(ofn,"w")

outputfile.write("#TumourCov\t"+str(C_cov)+"\n")
outputfile.write("#NormalCov\t"+str(N_cov)+"\n")
outputfile.write("#AberrantCellFraction\t"+str(acf)+"\n")
outputfile.write("#Ploidy\t"+str(ploidy)+"\n")

while line1:
		line1_split=line1.rstrip().split("\t")
		line2_split=line2.rstrip().split("\t")

		ch1=(line1_split[0])
		ch2=(line2_split[0])

		pos1=int(line1_split[1])
		try:
			pos2=int(line2_split[1])
		except:
			print(line1)
			print(line2)
			sys.exit()

		if ch1!=ch2 or pos1!=pos2:
				print "different position!"
				print line1
				print line2
				if ch1 < ch2:
					line1=inputfile1.readline()
					continue
				elif ch2 < ch1:
					line2=inputfile2.readline()
					continue
				elif pos1<pos2:
					line1=inputfile1.readline()
					continue
				elif pos2<pos1:
					line2=inputfile2.readline()
					continue
				else:
					break
		try:
			cov1=float(line1_split[-1])
			cov2=float(line2_split[-1])
		except:
			outputfile.write(line1_split[0]+"\t"+str(pos1)+"\t"+str(line1_split[-1])+"\t"+str(line2_split[-1])+"\t"+"NA\n")
			line1=inputfile1.readline()
			line2=inputfile2.readline()
			continue
		abs_CN=absolute_copy(cov1,cov2)

		if abs_CN >0:
				outputfile.write(line1_split[0]+"\t"+str(pos1)+"\t"+str(cov1)+"\t"+str(cov2)+"\t"+str(abs_CN)+"\n")
		else:
				outputfile.write(line1_split[0]+"\t"+str(pos1)+"\t"+str(cov1)+"\t"+str(cov2)+"\t"+"NA\n")

		line1=inputfile1.readline()
		line2=inputfile2.readline()
print "done"


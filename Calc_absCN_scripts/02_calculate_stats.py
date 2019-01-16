import sys
fn=sys.argv[1]
inputfile=file(fn)
line=inputfile.readline()
while "#"==line[0]:
	line=inputfile.readline()
line=inputfile.readline() #header

ofn=fn+".covstat"
outputfile=file(ofn,"w")
outputfile.write("filename\tregion\tthroughput\tavg_depth\n")

covered_region=0
covered_depth=0
avg_depth=0
while line:
	line_split=line.rstrip().split("\t")
#	print(line_split)
	this_covered_region=int(line_split[3])
	this_covered_depth=int(line_split[5])
	covered_region+=this_covered_region
	covered_depth+=this_covered_depth

	line=inputfile.readline()

avg_depth=round(covered_depth/float(covered_region),2)

outputfile.write(fn+"\t"+str(covered_region)+"\t"+str(covered_depth)+"\t"+str(avg_depth)+"\n")

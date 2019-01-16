# Lung Adenocarcinoma Structural Variation and Amplification Timing Analyses

---
##1. Prepare BAMs
We prepared paired tumor BAMs and normal BAMs. (e.g. A.tumor.bam, A.normal.bam). BAM files were aligned to the human genome reference (GRCh37 without the “chr” prefix in the contig name) by BWA-mem.

##2. Run Delly somatic call
We ran Delly v0.7.6 (structural variation caller) with command as below. 
Usage:
	delly call –t DEL –q 15 –n –o <Output file name> -g <Reference fasta> <Tumor BAM> <Normal BAM>
	delly call –t <DUP, INV, or TRA> -q 15 –o <Output file name> -g <Reference fasta> <Tumor BAM> <Normal BAM>
Example of command line:
	delly call –t DEL –q 15 –n –o test.DEL.bcf –g ref.fa test.tumor.bam test.normal.bam
	delly call –t DUP –q 15 –o test.DUP.bcf –g ref.fa test.tumor.bam test.normal.bam
	delly call –t INV –q 15 –o test.INV.bcf –g ref.fa test.tumor.bam test.normal.bam
	delly call –t TRA –q 15 –o test.TRA.bcf –g ref.fa test.tumor.bam test.normal.bam
Example of output file name:
	test.DEL.bcf, test.DUP.bcf, test.INV.bcf, test.TRA.bcf

##3. Merge outputs
To allow the next steps to be easily implemented, we merged multiple BCF output files into a VCF file for each sample.
Usage:
	bcftools concat –a –O v –o <Output.vcf> <Delly BCFs (output of step 2)>
Example of command line:
	bcftools concat –a –O v –o test.delly.vcf test.DEL.bcf test.DUP.bcf test.INV.bcf test.TRA.bcf
Example of output file name:
	test.delly.vcf

##4. Make panel of normal
We made a panel of normal (PON) file by merging multiple Delly VCFs. 
Usage:
	python Making_PON_Delly.py <Input text>
Format of Input text:
	Tab-delimited Delly VCF(output of step 3) list
	```
		ID1	/path/to/A.vcf
		ID2	/path/to/B.vcf
		ID3	/path/to/C.vcf
		…
	```
Example of command line:
	python Making_PON_Delly.py Delly_VCFs_list.txt
Example of output file name:
	PON.delly.txt

##5. SV processing and annotation
To distinguish true positives from false positives in the next filtering step, we used in-house scripts to annotate multiple columns to the Delly VCFs. These processes were done by a shell script running a series of in-house Python scripts. The individual steps of the process are described below.
Usage:
	sh Delly_annotation.sh <Input Delly VCF (output of step 3)> <column number of tumor in Input Delly VCF (10 or 11)> <Tumor BAM> <Normal BAM> <PON file> <DIR of SV_annot_scripts>
Example of command line:
	sh Delly_annotation.sh test.delly.vcf 10 test.tumor.bam test.normal.bam PON.delly.txt /path/to/Delly_annotation_scripts
Example of output file name:
	test.delly.vcf.somatic.annotated

##6. Filter SVs
Usage:
	python annotated_SV_filter.py <Input SV file (output of step 5)>
Example of command line:
	python annotated_SV_filter.py test.delly.vcf.somatic.annotated
Example of output file name:
	test.delly.vcf.somatic.annotated.fi

##7. Add breakpoints and edit columns
Usage:
	python annotated_SV_BPadd_edit.py <Input SV file (output of step 6)> <Normal BAM>
Example of command line:
	python annotated_SV_BPadd_edit.py test.delly.vcf.somatic.annotated.fi test.normal.bam
Example of output file name:
	test.delly.vcf.somatic.anotated.fi.BPedit
##8. Check all SVs using integrative genome browser (IGV)
##9. Clustering SVs
Usage:
	sh SV_clustering.sh <Input SV file (output of step 7 or 8)> <DIR of SV_cluster_scripts)>
Example of command line:
	sh SV_clustering.sh test.delly.vcf.somatic.annotated.fi.BPedit /path/to/SV_clustering_scripts
Example of output file name:
	test.delly.vcf.somatic.annotated.fi.BPedit.clustered

##10. Classification of complex clusters
10-1. Calculate absCN of each 100Kb bin
Usage:
	sh Get_100kb_absCN.sh <Tumor pileup file> <Normal pileup file> <Cellularity> <Ploidy> <gender (XX or XY)> <DIR of Calc_absCN_scripts>
Example of command line:
	sh Get_100kb_absCN.sh test.tumor.pileup.gz test.normal.pileup.gz 0.62 3.6 XY /path/to/Calc_absCN_scripts
Example of output file:
	test.tumor.pileup.100kbcov.absCN.gen_fi
	test.tumor.pileup.100kbcov.absCN.gen_fi.chrCN
10-2. Classify complex clusters
Usage:
	sh Classify_complexSV.sh <input SV file (output of step 9)> <100kb bin AbsCN file (output of step 10-1)> <chrCN file (output of step 10-1)> <DIR of scripts> <reference fasta index file>
Example of command line:
	sh Classify_complexSV.sh test.delly.vcf.somatic.annotated.fi.BPedit.clustered test.tumor.pileup.100kbcov.absCN.gen_fi test.tumor.pileup.100kbcov.absCN.gen_fi.chrCN /path/to/Classify_complexSV_scripts ref.fa.fai
Example of output file name:
	test.delly.vcf.somatic.annotated.fi.BPedit.clustered.gap_seg.100kbAbsCN.complex_class

##11. Amplification timing analysis
11-1. Merge SV breakpoints with segments of CNV
Usage:
	sh Merge_SV_CNV.sh <CNV segments.txt (output of Sequenza)> <Input SV file (output of step 8)> <DIR of scripts>
Example of command line:
	sh Merge_SV_CNV.sh test.segments.txt test.delly.vcf.somatic.annotated.fi.BPedit
Example of output file name:
	test.segments.txt.clean.SV_CNV_bp.txt

11-2. Calculate expected number of early SNVs in amplified segments
	Requirements:
		Input SNV file should have columns as below.
			1st: chromosome without “chr” prefix
			2nd: position
			9th: reference read count 
			10th: altered read count
	Usage:
		sh Calc_early_SNV.sh <CNV segments.txt (output of Sequenza)> <Input SNV file> <Input SV file (output of step 9)> <purity> <DIR of scripts> <sample id> <DIR of output>

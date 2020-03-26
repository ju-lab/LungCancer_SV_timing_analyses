#!/bin/bash
vcf=$1
tumor_col=$2
tumorBam=$3
normalBam=$4
pon=$5
srcDir=$6
outDir=$(dirname $vcf)


# START
echo $vcf > $log
echo "Starting:filter somatic"
(python $srcDir/01.filter_somatic_delly.py $vcf $tumor_col) || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:sorting"
(python $srcDir/02.sorting_delly.py $vcf.somatic) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic
echo "done"
echo "Starting:annotate PON"
(python $srcDir/03.annotate_PON.py $vcf.somatic.sort $pon $srcDir/hg19_ref.fa.fai) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic.sort
echo "done"
echo "Starting:find BP"
(python $srcDir/04.find_BP.py $vcf.somatic.sort.pon $tumorBam $normalBam) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic.sort.pon
echo "done"
echo "Starting:BP adjustment"
(python $srcDir/05.BP_adjustment.py $vcf.somatic.sort.pon.BPinfo) || { c=$?;echo "Error";exit $c; } 
rm $vcf.somatic.sort.pon.BPinfo
echo "done"
echo "Starting:count fragments and find newBP"
(python $srcDir/06.07.count_frag_find_newBP.py $vcf.somatic.sort.pon.BPinfo.BPadj $tumorBam $normalBam) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic.sort.pon.BPinfo.BPadj
echo "done"
echo "Starting:annotate background information"
(python $srcDir/08.annotate_BGinfo.py $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf $tumorBam $normalBam) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf
echo "done"
echo "Starting:annotate paired normal same clipping"
(python $srcDir/09.PN_same_clipping.py $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo $tumorBam $normalBam) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo
echo "done"
echo "Starting:annotate mapq starting position variance"
(python $srcDir/10.MAPQ_start_pos.py $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo.pnsc $tumorBam) || { c=$?;echo "Error";exit $c; }
rm $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo.pnsc
mv $vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo.pnsc.mqpos $vcf.somatic.annotated
echo "All Done"

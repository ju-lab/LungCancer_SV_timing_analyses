#!/bin/bash
input=$1
absCN=$2
chrCN=$3
srcDir=$4
reffai=$5
outDir=$(dirname $input)

log=$outDir/$input.Classify.log

# START
echo $input > $log
echo "Starting:Gap Seg annot"
(python $srcDir/01_SV_gap_seg_annot.py $input $reffai) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:100kbAbsCN annot"
(python $srcDir/02_SV_100kb_absCN_annot.py $input.gap_seg $absCN ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $input.gap_seg
echo "Starting:Classification"
(python $srcDir/03_complexSV_classify.py $input.gap_seg.100kbAbsCN $chrCN) &>> $log || { c=$?;echo "Error";exit $c; }
rm $input.gap_seg.100kbAbsCN
echo "All Done"

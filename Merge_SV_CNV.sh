#!/bin/bash
Seqz=$1
SVfn=$2
srcDir=$3
outDir=$(dirname $Seqz)

log=$outDir/$Seqz.Classify.log

# START
echo $input > $log
echo "Starting:Sequenza edit"
(python $srcDir/01_seqz_segment_clear.py $Seqz) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Make merged breakpoints"
(python $srcDir/02_merge_SV_CNV_BP.py $Seqz.clean $SVfn ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "All done"

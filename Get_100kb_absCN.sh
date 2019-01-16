#!/bin/bash
tpileup=$1
npileup=$2
cellularity=$3
ploidy=$4
gender=$5
srcDir=$6
outDir=$(dirname $tpileup)

log=$outDir/$tpileup.100kbAbsCN.log

# START
echo $input > $log
echo "Starting:Tumor binning"
(python $srcDir/01_get_coverage.py $tpileup) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Normal binning"
(python $srcDir/01_get_coverage.py $npileup) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Calculate mean depth"
tpileup=$(echo $tpileup | sed 's/.gz//')
npileup=$(echo $npileup | sed 's/.gz//')
(python $srcDir/02_calculate_stats.py $tpileup.100kbcov) &>> $log || { c=$?;echo "Error";exit $c; }
(python $srcDir/02_calculate_stats.py $npileup.100kbcov) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Calculate absolute CN"
mTdp=$(sed '1d;s/.*\t//' $tpileup.100kbcov.covstat)
mNdp=$(sed '1d;s/.*\t//' $npileup.100kbcov.covstat)
(python $srcDir/03_report_Tspecific_absolute_CN.py $tpileup.100kbcov $npileup.100kbcov $mTdp $mNdp $cellularity $ploidy ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Adjust X Y CN"
(python $srcDir/04_absCN_gender_adj.py $tpileup.100kbcov.absCN $gender) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Calculate median chromosomal CN"
(python $srcDir/05_calc_median_chr_CN.py $tpileup.100kbcov.absCN.gen_fi) &>> $log || { c=$?;echo "Error";exit $c; }
rm $tpileup.100kbcov
rm $npileup.100kbcov
rm $tpileup.100kbcov.covstat
rm $npileup.100kbcov.covstat
rm $tpileup.100kbcov.absCN
echo "All Done"


#!/bin/bash
input=$1
srcDir=$2
outDir=$(dirname $input)

log=$outDir/$input.SVcluster.log

# START
echo $input > $log
echo "Starting:clustering"
(python $srcDir/01.SV_cluster_proximity.py $input) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:split_simple"
(python $srcDir/02.split_simple_from_cluster.py $input.cluster) &>> $log || { c=$?;echo "Error";exit $c; }
rm $input.cluster
echo "done"
echo "Starting:edit table"
(python $srcDir/03.edit_to_table.py $input.cluster.split_simple) &>> $log || { c=$?;echo "Error";exit $c; }
rm $input.cluster.split_simple
mv $input.cluster.split_simple.edit $input.clustered
echo "All Done"

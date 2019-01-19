#!/bin/bash
Seqz=$1
SNVfn=$2
SVfn=$3
purity=$4
srcDir=$5
sampleid=$6
outDir=$7

log=$outDir/$sampleid.calc_early.log

# START
echo $input > $log
echo "Starting:Sequenza edit"
(python $srcDir/01_seqz_segments_edit.py $Seqz) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:mark broad amplification"
(python $srcDir/02_detect_broad_amp.py $Seqz.edit $srcDir/hg19_cytoBand.txt ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $Seqz.edit
echo "Starting:Annotate CNV to SNV"
(python $srcDir/03_SNV_seqz_edit_annot.py $SNVfn $Seqz.edit.broadAMP ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Annotate mutCN to SNV"
(python $srcDir/04_SNV_mutCN_CCF_annot.py $SNVfn.cninfo $purity ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $SNVfn.cninfo
echo "Starting:SNV file sort"
(python $srcDir/04.2_vcf_sorting.py $SNVfn.cninfo.scF) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $SNVfn.cninfo.scF
echo "Starting:Annotate kataegis to SNV"
(python $srcDir/05_SNV_annot_kataegis.py $SNVfn.cninfo.scF.sorted) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $SNVfn.cninfo.scF.sorted
echo "Starting:Annotate SNV probability"
(python $srcDir/06_SNV_timing_probability.py $SNVfn.cninfo.scF.sorted.kat $purity) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
rm $SNVfn.cninfo.scF.sorted.kat
echo "Starting:Annotate earlySNV to CNV"
(python $srcDir/07_Annot_EarlySNVCount_to_AmpSeg.py $Seqz.edit.broadAMP $SNVfn.cninfo.scF.sorted.kat.proba $sampleid $outDir) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Annotate earlySNV of cluster"
(python $srcDir/08_cluster_earlySNV_annot.py $outDir/$sampleid.ampseg.earlysnv.txt $SVfn $sampleid $outDir ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Merge earlySNV of cluster"
(python $srcDir/09_cluster_earlySNV_merge.py $outDir/$sampleid.cl_earlysnv ) &>> $log || { c=$?;echo "Error";exit $c; }
echo "done"
echo "Starting:Annotate earlySNV of broad amplification"
(python $srcDir/10_broad_earlySNV_annot.py $outDir/$sampleid.ampseg.earlysnv.txt $sampleid $outDir) &>> $log || { c=$?;echo "Error";exit $c; }
rm $outDir/$sampleid.ampseg.earlysnv.txt
rm $outDir/$sampleid.cl_earlysnv
echo "All done"

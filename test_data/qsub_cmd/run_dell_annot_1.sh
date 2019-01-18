#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/sypark/github_scripts/LungCancer_SV_timing/test_data/qsub_sdout/run_dell_annot_1.sh.sdout
cd /home/users/sypark/github_scripts/LungCancer_SV_timing/test_data
sh ../Delly_annotation.sh test.delly.vcf 10 ~/00_Project/05_Rearrangement/02_Bam/LU-F16.tumor.bam ~/00_Project/05_Rearrangement/02_Bam/LU-F16.normal.bam ~/00_Project/05_Rearrangement/19_revision_total_samples/07_delly/01_normal_dic/merged_144s_delly_normal_panel.txt ../Delly_annotation_scripts/ &> run_delly.log
wait

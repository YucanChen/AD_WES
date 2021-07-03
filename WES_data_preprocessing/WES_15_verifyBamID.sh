#!/bin/bash

cd /WES_output_combine/
for i in /recall/*_dedup_recal.bam
do
j=${i##*/WES_output/recall/}
verifyBamID --vcf ./filtered_vcf/WES_filter_pass_new.snp.vcf --bam $i --out /WES_output_combine/verifyBamID_files/${j%_dedup_recal.bam}_verifyBamID --verbose --ignoreRG &
done

for i in /recall/*_dedup_recal.bam
do
j=${i##*/WES_output_batch2/recall/}
verifyBamID --vcf /filtered_vcf/WES_filter_pass_new.snp.vcf --bam $i --out /WES_output_combine/verifyBamID_files/${j%_dedup_recal.bam}_verifyBamID --verbose --ignoreRG &
done
wait

exit 0

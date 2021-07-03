#!/bin/bash

for i in /recall/*_dedup_recal.bam
do
j=${i##*/recall/}
samtools flagstat -@ 16 $i > /WES_output_combine/afterrecall_QC/samtools_flagstat/${j%_dedup_recal.bam}_flagstat.txt &
done

for i in /recall/*_dedup_recal.bam
do
j=${i##*/recall/}
samtools flagstat -@ 16 $i > /WES_output_combine/afterrecall_QC/samtools_flagstat/${j%_dedup_recal.bam}_flagstat.txt &
done
wait

exit 0

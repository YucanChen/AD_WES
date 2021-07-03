#!/bin/bash
for i in /bam/*.bam
do
(samtools flagstat -@ 16 $i > /bam/bam_stat/$(basename $i _bwa.bam).bam_stat) &
done
wait
exit 0

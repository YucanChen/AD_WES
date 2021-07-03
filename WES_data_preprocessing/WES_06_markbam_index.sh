#!/bin/bash
cd /bam/
for i in /bam/*_marked.bam
do
(samtools index $i) &
done
wait
exit 0

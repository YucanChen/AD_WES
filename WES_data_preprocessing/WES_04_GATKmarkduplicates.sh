#!/bin/bash
for mapFile in /bam/*_bwa.bam
do
(~/gatk-4.1.0.0/gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" MarkDuplicates -I ${mapFile} -O ${mapFile%_bwa.bam}_marked.bam -M ${mapFile%_bwa.bam}.duplicate_metrics) &
done
wait
exit 0

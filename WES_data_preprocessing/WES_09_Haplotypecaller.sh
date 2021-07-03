#!/bin/bash
GATK=~/gatk-4.1.0.0/gatk
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
hg38exon_interval=~/gatk-4.1.0.0/bundle_hg38/targets.preprocessed.interval.list
for i in /recall/*_dedup_recal.bam
do
j=${i##*/recall/}
outfile=/gvcf/${j%_dedup_recal.bam}_dedup_recal.g.vcf
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller -R $ref -I $i -pairHMM AVX_LOGLESS_CACHING -L $hg38exon_interval --native-pair-hmm-threads 8 -O $outfile -ERC GVCF -stand-call-conf 10 &
done
wait

exit 0

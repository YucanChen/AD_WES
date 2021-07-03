#!/bin/bash
for i in /cleandata/R19*_1P.fastq.gz
do
filename=$(basename $i _1P.fastq.gz)
ref=~/gatk-4.1.0.0/bundle_hg38/gatk_hg38
mapFile=/bam/${filename}_bwa.bam
bwa mem -t 16 -M -R @RG\\tID:${filename}\\tSM:${filename}\\tPL:illumina\\tLB:${filename} $ref $i ${i%1P.fastq.gz}2P.fastq.gz |samtools sort - -@ 16 -m 4G -T ${i%_1P.fastq.gz} -o $mapFile
done

for i in /cleandata/R2000*_1P.fastq.gz
do
filename=$(basename $i _1P.fastq.gz)
ref=~/gatk-4.1.0.0/bundle_hg38/gatk_hg38
mapFile=/bam/${filename}_bwa.bam
bwa mem -t 16 -M -R @RG\\tID:${filename}\\tSM:${filename}\\tPL:illumina\\tLB:${filename} $ref $i ${i%1P.fastq.gz}2P.fastq.gz |samtools sort - -@ 16 -m 4G -T ${i%_1P.fastq.gz} -o $mapFile
done
exit 0

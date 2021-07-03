#!/bin/bash
snp=~/gatk-4.1.0.0/bundle_hg38/dbsnp_146.hg38.vcf.gz
indel=~/gatk-4.1.0.0/bundle_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
snp1000G=~/gatk-4.1.0.0/bundle_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
hg38exon_interval=~/gatk-4.1.0.0/bundle_hg38/targets.preprocessed.interval.list
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
GATK=~/gatk-4.1.0.0/gatk

for mapFile in /bam/*_marked.bam
do
j=${mapFile##*/bam/}
outfile=/recall/${j%_marked.bam}_recal.table
$GATK --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator -L $hg38exon_interval -R $ref -I $mapFile --known-sites $indel --known-sites $snp1000G --known-sites $snp -O $outfile &
done
wait

exit 0

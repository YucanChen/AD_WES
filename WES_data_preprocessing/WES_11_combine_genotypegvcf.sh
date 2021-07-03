#!/bin/bash

####CombineGVCFs
GATK=~/gatk-4.1.0.0/gatk
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
output=/combine_genotype_vcf/WES_combine.g.vcf
$GATK CombineGVCFs -R $ref $(for i in /WES_output_combine/gvcf_filtered/*_dedup_recal_highDPrm.g.vcf;do echo "--variant $i";done) -O $output

####genotypeGVCFs
# genotype gvcf files
hg38exon_interval=~/gatk-4.1.0.0/bundle_hg38/targets.preprocessed.interval.list
infile=/WES_output_combine/combine_genotype_vcf/WES_combine.g.vcf
outfile=/WES_output_combine/combine_genotype_vcf/WES_combine.vcf
$GATK --java-options "-Xmx40G -XX:+UseParallelGC -XX:ParallelGCThreads=10" GenotypeGVCFs -R $ref -V $infile -L $hg38exon_interval -O $outfile

exit 0

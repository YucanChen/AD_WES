#!/bin/bash

GATK=~/gatk-4.1.0.0/gatk
vcfFile_recal_v2=/WES_output_combine/VQSR/WES_combine.recal_v2.vcf
bcftools_tmp1=/WES_output_combine/VQSR/WES_combine.recal_v2_bcf1.vcf
bcftools_tmp2=/WES_output_combine/VQSR/WES_combine.recal_v2_bcf2.vcf
indel_filtered=/WES_output_combine/VQSR/WES_distance_filter.indel.vcf
snp_filtered=/WES_output_combine/VQSR/WES_distance_filter.snp.vcf
snp_tmp=/WES_output_combine/VQSR/WES_distance_tmp.snp.vcf
indel_tmp=/WES_output_combine/VQSR/WES_distance_tmp.indel.vcf
filtered_indel_pass=/WES_output_combine/filtered_vcf/WES_filter_pass_new.indel.vcf
filtered_snp_pass=/WES_output_combine/filtered_vcf/WES_filter_pass_new.snp.vcf

####根据距离筛选indels和snps
bcftools filter --IndelGap 50 $vcfFile_recal_v2 --output-type z --output $bcftools_tmp1.gz
bcftools filter --SnpGap 100 $bcftools_tmp1.gz --output-type z --output $bcftools_tmp2.gz
gunzip $bcftools_tmp2.gz
# 使用SelectVariants，选出SNP
$GATK SelectVariants -select-type SNP -V $bcftools_tmp2 -O $snp_filtered
# 使用SelectVariants，选出Indel
$GATK SelectVariants -select-type INDEL -V $bcftools_tmp2 -O $indel_filtered
# Filtering
$GATK VariantFiltration -V $snp_filtered --filter-expression "QD<2.0 || InbreedingCoeff<-0.2 || AC=0" --filter-name "Filter" -O $snp_tmp
$GATK VariantFiltration -V $indel_filtered --filter-expression "QD<2.0 || InbreedingCoeff<-0.2 || AC=0" --filter-name "Filter" -O $indel_tmp

$GATK SelectVariants --exclude-filtered true -V $snp_tmp -O $filtered_snp_pass
$GATK SelectVariants --exclude-filtered true -V $indel_tmp -O $filtered_indel_pass

exit 0

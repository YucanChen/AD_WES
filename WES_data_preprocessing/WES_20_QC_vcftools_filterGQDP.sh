#!/bin/bash

snps_vcftools=/WES_output_combine/QC_passed_vcf/WES_qc_pass.snp.vcf
indels_vcftools=/WES_output_combine/QC_passed_vcf/WES_qc_pass.indel.vcf
snp_GQDP_2bfilter=/WES_output_combine/QC_passed_vcf/vcftools_merge_snp_GQDP_2bfilter.txt
indel_GQDP_2bfilter=/WES_output_combine/QC_passed_vcf/vcftools_merge_indel_GQDP_2bfilter.txt
snp_qc_final=/WES_output_combine/QC_passed_vcf/WES_qc_pass_final.snp.vcf
indel_qc_final=/WES_output_combine/QC_passed_vcf/WES_qc_pass_final.indel.vcf

vcftools --vcf $snps_vcftools --exclude-positions $snp_GQDP_2bfilter --stdout --recode --recode-INFO-all > $snp_qc_final &
vcftools --vcf $indels_vcftools --exclude-positions $indel_GQDP_2bfilter --stdout --recode --recode-INFO-all > $indel_qc_final &
wait

exit 0

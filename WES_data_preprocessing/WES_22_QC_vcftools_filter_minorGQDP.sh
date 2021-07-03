#!/bin/bash

snps_vcftools=/WES_output_combine/QC_passed_vcf/WES_qc_pass_final.snp.vcf
indels_vcftools=/WES_output_combine/QC_passed_vcf/WES_qc_pass_final.indel.vcf
snp_2bleft_minor=/WES_output_combine/QC_passed_vcf/snp_2bleft_minor.txt
indel_2bleft_minor=/WES_output_combine/QC_passed_vcf/indel_2bleft_minor.txt
snp_final_pass=/WES_output_combine/QC_passed_vcf/WES_qc_filtered.snp.vcf
indel_final_pass=/WES_output_combine/QC_passed_vcf/WES_qc_filtered.indel.vcf

vcftools --vcf $snps_vcftools --positions $snp_2bleft_minor --stdout --recode --recode-INFO-all > $snp_final_pass &
vcftools --vcf $indels_vcftools --positions $indel_2bleft_minor --stdout --recode --recode-INFO-all > $indel_final_pass &
wait

exit 0

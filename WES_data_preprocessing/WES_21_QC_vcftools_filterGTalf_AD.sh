#!/bin/bash

snp_qc_final=/WES_output_combine/QC_passed_vcf/WES_qc_pass_final.snp.vcf
indel_qc_final=/WES_output_combine/QC_passed_vcf/WES_qc_pass_final.indel.vcf
snp_prefix=/WES_output_combine/QC_passed_vcf/snp_GQDP_filtered.out
indel_prefix=/WES_output_combine/QC_passed_vcf/indel_GQDP_filtered.out

vcftools --extract-FORMAT-info GT --vcf $snp_qc_final --out $snp_prefix &
vcftools --extract-FORMAT-info GT --vcf $indel_qc_final --out $indel_prefix &
vcftools --extract-FORMAT-info AD --vcf $snp_qc_final --out $snp_prefix &
vcftools --extract-FORMAT-info AD --vcf $indel_qc_final --out $indel_prefix &
vcftools --extract-FORMAT-info GQ --vcf $snp_qc_final --out $snp_prefix &
vcftools --extract-FORMAT-info GQ --vcf $indel_qc_final --out $indel_prefix &
wait

exit 0

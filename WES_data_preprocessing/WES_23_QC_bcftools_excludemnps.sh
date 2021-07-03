#!/bin/bash

snp_final_pass=/WES_output_combine/QC_passed_vcf/WES_qc_filtered.snp.vcf
indel_final_pass=/WES_output_combine/QC_passed_vcf/WES_qc_filtered.indel.vcf
seletype_snp=/WES_output_combine/QC_passed_vcf/WES_type_selected.snp.vcf
seletype_indel=/WES_output_combine/QC_passed_vcf/WES_type_selected.indel.vcf
bcftools view --exclude-types mnps,ref $snp_final_pass -Ov > $seletype_snp &
bcftools view --exclude-types mnps,ref $indel_final_pass -Ov > $seletype_indel &
wait

exit 0

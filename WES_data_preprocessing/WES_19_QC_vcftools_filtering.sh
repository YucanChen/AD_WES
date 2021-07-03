#!/bin/bash

snps=/WES_output_combine/filtered_vcf/WES_filter_pass_new.snp.vcf
indels=/WES_output_combine/filtered_vcf/WES_filter_pass_new.indel.vcf
masked_bed=/bed/hg38.masked.bed
rmmask_snps=/WES_output_combine/filtered_vcf/WES_rmmasked.snp.vcf
rmmask_indels=/WES_output_combine/filtered_vcf/WES_rmmasked.indel.vcf
tmp_snps=/WES_output_combine/QC_passed_vcf/WES_tmp.snp.vcf
tmp_indels=/WES_output_combine/QC_passed_vcf/WES_tmp.indel.vcf
snps_vcftools=/WES_output_combine/QC_passed_vcf/WES_qc_pass.snp.vcf
indels_vcftools=/WES_output_combine/QC_passed_vcf/WES_qc_pass.indel.vcf
remove_list=/WES_output_combine/afterrecall_QC/fitered_individuals_relatandpca.txt
snp_prefix=/WES_output_combine/QC_passed_vcf/snp.out
indel_prefix=/WES_output_combine/QC_passed_vcf/indel.out

vcftools --vcf $snps --exclude-bed $masked_bed --stdout --recode --recode-INFO-all > $rmmask_snps &
vcftools --vcf $indels --exclude-bed $masked_bed --stdout --recode --recode-INFO-all > $rmmask_indels &
wait
#vcftools --vcf $rmmask_snps --remove $remove_list --stdout --recode --recode-INFO-all > $tmp_snps &
#vcftools --vcf $rmmask_indels --remove $remove_list --stdout --recode --recode-INFO-all > $tmp_indels &
wait
vcftools --vcf $rmmask_snps  --max-missing 0.95 --stdout --recode --recode-INFO-all > $snps_vcftools &
vcftools --vcf $rmmask_indels  --max-missing 0.95 --stdout --recode --recode-INFO-all > $indels_vcftools &
wait
vcftools --vcf $snps_vcftools --missing-indv --out $snp_prefix &
vcftools --vcf $snps_vcftools --missing-indv --out $indel_prefix &

vcftools --vcf $snps_vcftools --extract-FORMAT-info GQ --out $snp_prefix &
vcftools --vcf $indels_vcftools --extract-FORMAT-info GQ --out $indel_prefix &
vcftools --vcf $snps_vcftools --extract-FORMAT-info DP --out $snp_prefix &
vcftools --vcf $indels_vcftools --extract-FORMAT-info DP --out $indel_prefix &
wait

exit 0

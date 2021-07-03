#!/bin/bash

input_vcf=/WES_output_combine/filtered_vcf/WES_filter_pass_new.snp.vcf
common_snps=/WES_output_combine/afterrecall_QC/PCA/WES_filter_pass_common05.snp.vcf
plinkout=/WES_output_combine/afterrecall_QC/PCA/WES_filter_pass_common05_pca
transfer_conf=/WES_output_combine/afterrecall_QC/PCA/transfer.conf
runningpca_conf=/WES_output_combine/afterrecall_QC/PCA/runningpca.conf
myfile_pca=/WES_output_combine/afterrecall_QC/PCA/WES_filter_pass_common05_plinkpca

# select common snps for PCA analyis
vcftools --maf 0.05 --vcf $input_vcf --recode --recode-INFO-all --stdout > $common_snps
# data format conversion
vcftools  --vcf $common_snps --plink --out $plinkout
# plink pca file preparation
plink --vcf $common_snps --make-bed --recode --out $plinkout --double-id
# plink pca calculation
plink --bfile $plinkout --pca 10 --out $myfile_pca

exit 0

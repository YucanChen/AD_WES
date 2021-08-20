#!/bin/bash

# merge indel and snp files seperately
# sort snps and indels seperately：
in_snp=~/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38.snp.snpEff.vcf
in_indel=~/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38.indel.snpEff.vcf
out_snp=~/burden_test/case_batch2_dementia/WES_highquality_hg38_sorted.snp.snpEff.vcf.gz
out_indel=~/burden_test/case_batch2_dementia/WES_highquality_hg38_sorted.indel.snpEff.vcf.gz
merged_variants=~/burden_test/case_batch2_dementia/WES_highquality_hg38_sorted.all.snpEff.vcf.gz
bcftools sort $in_snp -O z -o $out_snp && tabix -p vcf $out_snp &
bcftools sort $in_indel -O z -o $out_indel && tabix -p vcf $out_indel &
wait

# combine snps and indels：
bcftools concat $out_snp $out_indel -a -O z -o $merged_variants && tabix -p vcf $merged_variants

# Separating multi-allelics and left-aligning indels:
in_vcf=~/burden_test/case_batch2_dementia/WES_highquality_hg38_sorted.all.snpEff.vcf.gz
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
out_vcf=~/burden_test/case_batch2_dementia/WES_highquality_hg38_sorted.all.norm.vcf.gz
bcftools norm -m -any -f $ref $in_vcf | bgzip > $out_vcf && tabix -p vcf $out_vcf &
wait

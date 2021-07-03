#!/bin/bash

input_snp=/WES_output_combine/QC_passed_vcf/WES_type_selected.snp.vcf
input_indel=/WES_output_combine/QC_passed_vcf/WES_type_selected.indel.vcf
input_snp_bgz=/WES_output_combine/QC_passed_vcf/WES_type_selected.snp.vcf.bgz
input_indel_bgz=/WES_output_combine/QC_passed_vcf/WES_type_selected.indel.vcf.bgz
tmp_snpsubset_AD=/WES_output_combine/QC_passed_vcf/WES_type_selected_AD_tmp.snp.vcf.bgz
tmp_indelsubset_AD=/WES_output_combine/QC_passed_vcf/WES_type_selected_AD_tmp.indel.vcf.bgz
tmp_snpsubset_ADearly=/WES_output_combine/QC_passed_vcf/WES_type_selected_ADearly_tmp.snp.vcf.bgz
tmp_indelsubset_ADearly=/WES_output_combine/QC_passed_vcf/WES_type_selected_ADearly_tmp.indel.vcf.bgz
snpsubset_AD=/WES_output_combine/QC_passed_vcf/WES_type_selected_AD.snp.vcf
indelsubset_AD=/WES_output_combine/QC_passed_vcf/WES_type_selected_AD.indel.vcf
snpsubset_ADearly=/WES_output_combine/QC_passed_vcf/WES_type_selected_ADearly.snp.vcf
indelsubset_ADearly=/WES_output_combine/QC_passed_vcf/WES_type_selected_ADearly.indel.vcf

keep_list_AD=keep_list_AD.txt
keep_list_ADearly=keep_list_ADearly.txt

###Part0: bgzip tabix vcf files
cat $input_snp | bgzip -c > $input_snp_bgz && tabix -p vcf $input_snp_bgz &
cat $input_indel | bgzip -c > $input_indel_bgz && tabix -p vcf $input_indel_bgz
wait

###Part1: AD samples
# snp subset:
bcftools view -S $keep_list_AD $input_snp_bgz | bgzip -c > $tmp_snpsubset_AD && tabix -p vcf $tmp_snpsubset_AD &
# indel subset:
bcftools view -S $keep_list_AD $input_indel_bgz | bgzip -c > $tmp_indelsubset_AD && tabix -p vcf $tmp_indelsubset_AD &

###Part2: early AD samples
# snp subset:
bcftools view -S $keep_list_ADearly $input_snp_bgz | bgzip -c > $tmp_snpsubset_ADearly && tabix -p vcf $tmp_snpsubset_ADearly &
# indel subset:
bcftools view -S $keep_list_ADearly $input_indel_bgz | bgzip -c > $tmp_indelsubset_ADearly && tabix -p vcf $tmp_indelsubset_ADearly &
wait

###Part3: delete non-variant sites
bcftools view -c 1:nonmajor $tmp_snpsubset_AD > $snpsubset_AD &
bcftools view -c 1:nonmajor $tmp_indelsubset_AD > $indelsubset_AD &
bcftools view -c 1:nonmajor $tmp_snpsubset_ADearly > $snpsubset_ADearly &
bcftools view -c 1:nonmajor $tmp_indelsubset_ADearly > $indelsubset_ADearly &
wait

rm $tmp_snpsubset_AD && rm $tmp_indelsubset_AD && rm $tmp_snpsubset_ADearly && $tmp_indelsubset_ADearly && rm $tmp_snpsubset_AD.tbi && rm $tmp_indelsubset_AD.tbi && rm $tmp_snpsubset_ADearly.tbi && $tmp_indelsubset_ADearly.tbi

exit 0

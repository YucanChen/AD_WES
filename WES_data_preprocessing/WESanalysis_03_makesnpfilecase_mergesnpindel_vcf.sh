#!/bin/bash

# section 1: prepare makesnpfile
original_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.snpEff.vcf.gz
rare_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.rare.vcf.gz
PAV_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV.vcf.gz
cal_PAV_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D.vcf.gz
Syn_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.Syn.vcf.gz

original_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.snpEff_deXY.vcf.gz
rare_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.rare_deXY.vcf.gz
PAV_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_deXY.vcf.gz
cal_PAV_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.vcf.gz
Syn_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.Syn_deXY.vcf.gz

original_case_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.snpEff.vcf.gz
rare_case_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.rare.vcf.gz
PAV_case_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV.vcf.gz
Syn_case_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.Syn.vcf.gz

original_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.snpEff_deXY.vcf.gz
rare_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.rare_deXY.vcf.gz
PAV_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.vcf.gz
Syn_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.Syn_deXY.vcf.gz

java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $original_case "(ExAC_ALL<0.001) & (ExAC_AFR<0.001) & (ExAC_AMR<0.001) & (ExAC_EAS<0.001) & (ExAC_FIN<0.001) & (ExAC_NFE<0.001) & (ExAC_OTH<0.001) & (ExAC_SAS<0.001) & (esp6500siv2_ea<0.001) & (esp6500siv2_aa<0.001) & (esp6500siv2_all<0.001) & (AF_popmax<0.001) & (non_topmed_AF_popmax<0.001) & (non_neuro_AF_popmax<0.001) & (non_cancer_AF_popmax<0.001) & (controls_AF_popmax<0.001) & (X1000G_AFR<0.001) & (X1000G_ALL<0.001) & (X1000G_AMR<0.001) & (X1000G_EAS<0.001) & (X1000G_EUR<0.001) & (X1000G_SAS<0.001)" | bgzip -c > $rare_case && tabix -p vcf $rare_case && java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $rare_case "(ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'stop_gained')" | bgzip -c > $PAV_case && tabix -p vcf $PAV_case && java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $rare_case "ANN[*].EFFECT = 'synonymous_variant'" | bgzip -c > $Syn_case && tabix -p vcf $Syn_case
java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $PAV_case "(Polyphen2_HDIV_pred='D')" | bgzip -c > $cal_PAV_case && tabix -p vcf $cal_PAV_case

java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $original_case_indel "(ExAC_ALL<0.001) & (ExAC_AFR<0.001) & (ExAC_AMR<0.001) & (ExAC_EAS<0.001) & (ExAC_FIN<0.001) & (ExAC_NFE<0.001) & (ExAC_OTH<0.001) & (ExAC_SAS<0.001) & (esp6500siv2_ea<0.001) & (esp6500siv2_aa<0.001) & (esp6500siv2_all<0.001) & (AF_popmax<0.001) & (non_topmed_AF_popmax<0.001) & (non_neuro_AF_popmax<0.001) & (non_cancer_AF_popmax<0.001) & (controls_AF_popmax<0.001) & (X1000G_AFR<0.001) & (X1000G_ALL<0.001) & (X1000G_AMR<0.001) & (X1000G_EAS<0.001) & (X1000G_EUR<0.001) & (X1000G_SAS<0.001)" | bgzip -c > $rare_case_indel && tabix -p vcf $rare_case_indel && java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $rare_case_indel "(ANN[*].EFFECT has 'missense_variant') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'stop_gained')" | bgzip -c > $PAV_case_indel && tabix -p vcf $PAV_case_indel && java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $rare_case_indel "ANN[*].EFFECT = 'synonymous_variant'" | bgzip -c > $Syn_case_indel && tabix -p vcf $Syn_case_indel

vcftools --gzvcf $original_case --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $original_case_deXY && tabix -p vcf $original_case_deXY &
vcftools --gzvcf $rare_case --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $rare_case_deXY && tabix -p vcf $rare_case_deXY &
vcftools --gzvcf $PAV_case --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $PAV_case_deXY && tabix -p vcf $PAV_case_deXY &
vcftools --gzvcf $cal_PAV_case --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $cal_PAV_case_deXY && tabix -p vcf $cal_PAV_case_deXY &
vcftools --gzvcf $Syn_case --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $Syn_case_deXY && tabix -p vcf $Syn_case_deXY &

vcftools --gzvcf $original_case_indel --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $original_case_indel_deXY && tabix -p vcf $original_case_indel_deXY &
vcftools --gzvcf $rare_case_indel --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $rare_case_indel_deXY && tabix -p vcf $rare_case_indel_deXY &
vcftools --gzvcf $PAV_case_indel --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $PAV_case_indel_deXY && tabix -p vcf $PAV_case_indel_deXY &
vcftools --gzvcf $Syn_case_indel --not-chr chrX --not-chr chrY --stdout --recode --recode-INFO-all | bgzip -c > $Syn_case_indel_deXY && tabix -p vcf $Syn_case_indel_deXY &
wait

# section 2: merge snp and indel vcf
original_snp_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.vcf.gz
original_indel_case=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.vcf.gz
bgzip_snp_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.vcf.bgz
bgzip_indel_case=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.vcf.bgz

zcat $original_snp_case | bgzip -c > $bgzip_snp_case && tabix -p vcf $bgzip_snp_case &
zcat $original_indel_case | bgzip -c > $bgzip_indel_case && tabix -p vcf $bgzip_indel_case &
wait
merge_vcf_case=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY.merge.vcf.gz
bcftools concat -a $bgzip_snp_case $bgzip_indel_case -Oz -o $merge_vcf_case &
wait

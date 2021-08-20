#!/bin/bash

####case inner MAF adjust to：0.05
original_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.snpEff_deXY.vcf.gz
rare_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.rare_deXY.vcf.gz
Syn_case_indel_deXY=~/burden_test/case/WES_highquality_hg38_sorted.indel.Syn_deXY.vcf.gz

original_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.snpEff_deXY.vcf.gz
rare_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.rare_deXY.vcf.gz
PAV_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_deXY.vcf.gz
Syn_case_deXY=~/burden_test/case/WES_highquality_hg38_sorted.snp.Syn_deXY.vcf.gz

bgzip_snp_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.vcf.gz
bgzip_indel_case=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.vcf.gz
merge_vcf_case=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY.merge.vcf.gz

original_case_indel_deXY_tmp=~/burden_test/case/tmp1.vcf.gz
rare_case_indel_deXY_tmp=~/burden_test/case/tmp2.vcf.gz
Syn_case_indel_deXY_tmp=~/burden_test/case/tmp3.vcf.gz

original_case_deXY_tmp=~/burden_test/case/tmp4.vcf.gz
rare_case_deXY_tmp=~/burden_test/case/tmp5.vcf.gz
PAV_case_deXY_tmp=~/burden_test/case/tmp6.vcf.gz
Syn_case_deXY_tmp=~/burden_test/case/tmp7.vcf.gz

bgzip_snp_case_tmp=~/burden_test/case/tmp8.vcf.gz
bgzip_indel_case_tmp=~/burden_test/case/tmp9.vcf.gz
merge_vcf_case_tmp=~/burden_test/case/tmp10.vcf.gz

original_case_indel_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.snpEff_deXY_inmaf.vcf.bgz
rare_case_indel_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.rare_deXY_inmaf.vcf.bgz
Syn_case_indel_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.Syn_deXY_inmaf.vcf.bgz

original_case_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.snpEff_deXY_inmaf.vcf.bgz
rare_case_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.rare_deXY_inmaf.vcf.bgz
PAV_case_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_deXY_inmaf.vcf.bgz
Syn_case_deXY_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.Syn_deXY_inmaf.vcf.bgz

bgzip_snp_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.vcf.bgz
bgzip_indel_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.vcf.bgz
merge_vcf_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY_inmaf.merge.vcf.gz

homalt_chr11_2bfilter=~/burden_test/case/homalt_chr11_2bfilter.txt

vcftools --gzvcf $original_case_indel_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $original_case_indel_deXY_tmp &
vcftools --gzvcf $rare_case_indel_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $rare_case_indel_deXY_tmp &
vcftools --gzvcf $Syn_case_indel_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $Syn_case_indel_deXY_tmp &
vcftools --gzvcf $original_case_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $original_case_deXY_tmp &
vcftools --gzvcf $rare_case_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $rare_case_deXY_tmp &
vcftools --gzvcf $PAV_case_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $PAV_case_deXY_tmp &
vcftools --gzvcf $Syn_case_deXY --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $Syn_case_deXY_tmp &
vcftools --gzvcf $bgzip_snp_case --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $bgzip_snp_case_tmp &
vcftools --gzvcf $bgzip_indel_case --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $bgzip_indel_case_tmp &
vcftools --gzvcf $merge_vcf_case --exclude-positions $homalt_chr11_2bfilter --stdout --recode --recode-INFO-all | gzip > $merge_vcf_case_tmp &
wait

vcftools --gzvcf $original_case_indel_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $original_case_indel_deXY_inner && tabix -p vcf $original_case_indel_deXY_inner &
vcftools --gzvcf $rare_case_indel_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $rare_case_indel_deXY_inner && tabix -p vcf $rare_case_indel_deXY_inner &
vcftools --gzvcf $Syn_case_indel_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $Syn_case_indel_deXY_inner && tabix -p vcf $Syn_case_indel_deXY_inner &
vcftools --gzvcf $original_case_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $original_case_deXY_inner && tabix -p vcf $original_case_deXY_inner &
vcftools --gzvcf $rare_case_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $rare_case_deXY_inner && tabix -p vcf $rare_case_deXY_inner &
vcftools --gzvcf $PAV_case_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $PAV_case_deXY_inner && tabix -p vcf $PAV_case_deXY_inner &
vcftools --gzvcf $Syn_case_deXY_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $Syn_case_deXY_inner && tabix -p vcf $Syn_case_deXY_inner &
vcftools --gzvcf $bgzip_snp_case_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $bgzip_snp_case_inner && tabix -p vcf $bgzip_snp_case_inner &
vcftools --gzvcf $bgzip_indel_case_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $bgzip_indel_case_inner && tabix -p vcf $bgzip_indel_case_inner &
vcftools --gzvcf $merge_vcf_case_tmp --max-maf 0.05 --stdout --recode --recode-INFO-all | bgzip -c > $merge_vcf_case_inner && tabix -p vcf $merge_vcf_case_inner &
wait

rm $original_case_indel_deXY_tmp && rm $rare_case_indel_deXY_tmp && rm $Syn_case_indel_deXY_tmp && rm $original_case_deXY_tmp && rm $rare_case_deXY_tmp && rm $PAV_case_deXY_tmp && rm $Syn_case_deXY_tmp && rm $bgzip_snp_case_tmp && rm $bgzip_indel_case_tmp && rm $merge_vcf_case_tmp

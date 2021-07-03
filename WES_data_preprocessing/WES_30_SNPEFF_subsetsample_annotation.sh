#!/bin/bash

in_snp_AD=/WES_output_combine/ANNOVAR_annotation/WES_highquality_AD.snp.hg38_multianno.vcf
in_indel_AD=/WES_output_combine/ANNOVAR_annotation/WES_highquality_AD.indel.hg38_multianno.vcf
out_snp_AD=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38_AD.snp.snpEff.vcf
out_indel_AD=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38_AD.indel.snpEff.vcf

in_snp_ADearly=/WES_output_combine/ANNOVAR_annotation/WES_highquality_ADearly.snp.hg38_multianno.vcf
in_indel_ADearly=/WES_output_combine/ANNOVAR_annotation/WES_highquality_ADearly.indel.hg38_multianno.vcf
out_snp_ADearly=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38_ADearly.snp.snpEff.vcf
out_indel_ADearly=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38_ADearly.indel.snpEff.vcf

java -Xmx4g -jar ~/SnpEFF/snpEff/snpEff.jar GRCh38.86 $in_snp_AD > $out_snp_AD &
java -Xmx4g -jar ~/SnpEFF/snpEff/snpEff.jar GRCh38.86 $in_indel_AD > $out_indel_AD &
java -Xmx4g -jar ~/SnpEFF/snpEff/snpEff.jar GRCh38.86 $in_snp_ADearly > $out_snp_ADearly &
java -Xmx4g -jar ~/SnpEFF/snpEff/snpEff.jar GRCh38.86 $in_indel_ADearly > $out_indel_ADearly &
wait

exit 0

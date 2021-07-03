#!/bin/bash

in_snp=/WES_output_combine/ANNOVAR_annotation/WES_highquality.snp.hg38_multianno.vcf
in_indel=/WES_output_combine/ANNOVAR_annotation/WES_highquality.indel.hg38_multianno.vcf
out_snp=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38.snp.snpEff.vcf
out_indel=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38.indel.snpEff.vcf
java -Xmx4g -jar ~/SnpEFF/snpEff/snpEff.jar GRCh38.86 $in_snp > $out_snp &
java -Xmx4g -jar ~/SnpEFF/snpEff/snpEff.jar GRCh38.86 $in_indel > $out_indel &
wait

exit 0

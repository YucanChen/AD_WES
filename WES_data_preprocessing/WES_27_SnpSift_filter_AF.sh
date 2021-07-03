#!/bin/bash

in_snp=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38.snp.snpEff.vcf
in_indel=/WES_output_combine/SNPEFF_annotation/WES_highquality_hg38.indel.snpEff.vcf
filtered_snp=/WES_output_combine/SNPEFF_filtering/WES_highquality_hg38.snp.rare.vcf
filtered_indel=/WES_output_combine/SNPEFF_filtering/WES_highquality_hg38.indel.rare.vcf
java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $in_snp "(ExAC_ALL<0.001) & (ExAC_AFR<0.001) & (ExAC_AMR<0.001) & (ExAC_EAS<0.001) & (ExAC_FIN<0.001) & (ExAC_NFE<0.001) & (ExAC_OTH<0.001) & (ExAC_SAS<0.001) & (esp6500siv2_ea<0.001) & (esp6500siv2_aa<0.001) & (esp6500siv2_all<0.001) & (AF_popmax<0.001) & (non_topmed_AF_popmax<0.001) & (non_neuro_AF_popmax<0.001) & (non_cancer_AF_popmax<0.001) & (controls_AF_popmax<0.001) & (X1000G_AFR<0.001) & (X1000G_ALL<0.001) & (X1000G_AMR<0.001) & (X1000G_EAS<0.001) & (X1000G_EUR<0.001) & (X1000G_SAS<0.001)" > $filtered_snp &
java -Xmx4g -jar ~/SnpEFF/snpEff/SnpSift.jar filter -f $in_indel "(ExAC_ALL<0.001) & (ExAC_AFR<0.001) & (ExAC_AMR<0.001) & (ExAC_EAS<0.001) & (ExAC_FIN<0.001) & (ExAC_NFE<0.001) & (ExAC_OTH<0.001) & (ExAC_SAS<0.001) & (esp6500siv2_ea<0.001) & (esp6500siv2_aa<0.001) & (esp6500siv2_all<0.001) & (AF_popmax<0.001) & (non_topmed_AF_popmax<0.001) & (non_neuro_AF_popmax<0.001) & (non_cancer_AF_popmax<0.001) & (controls_AF_popmax<0.001) & (X1000G_AFR<0.001) & (X1000G_ALL<0.001) & (X1000G_AMR<0.001) & (X1000G_EAS<0.001) & (X1000G_EUR<0.001) & (X1000G_SAS<0.001)" > $filtered_indel &
wait

exit 0

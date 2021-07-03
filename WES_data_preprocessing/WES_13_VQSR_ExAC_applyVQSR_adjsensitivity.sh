#!/bin/bash

GATK=~/gatk-4.1.0.0/gatk
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
hg38exon_interval=~/gatk-4.1.0.0/bundle_hg38/targets.preprocessed.interval.list

vcfFile=/WES_output_combine/combine_genotype_vcf/WES_combine.vcf

output_recalSNP=/WES_output_combine/VQSR/WES_recalibrate_SNP.recal
output_tranchSNP=/WES_output_combine/VQSR/WES_recalibrate_SNP.tranches
output_rscriptSNP=/WES_output_combine/VQSR/output.SNP.recalibration_plots.rscript
recal_snpsindel=/WES_output_combine/VQSR/WES_recalibrated_snps_raw_indels.vcf

output_recalINDEL=/WES_output_combine/VQSR/WES_recalibrate_INDEL.recal
output_tranchINDEL=/WES_output_combine/VQSR/WES_recalibrate_INDEL.tranches
output_rscriptINDEL=/WES_output_combine/VQSR/output.INDEL.recalibration_plots.rscript

### the Ts/Tv obtained above is ~1.7, so I need to raise the sensitivity threshold for better practice:
vcfFile_recal_v2=/WES_output_combine/VQSR/WES_combine.recal_v2.vcf
# Apply recalibration to SNPs (SNP sensitivity 99.8)/Apply SNP VQSR on genotypes file
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=20" ApplyVQSR -V $vcfFile -R $ref -L $hg38exon_interval -O $recal_snpsindel  --recal-file $output_recalSNP --tranches-file $output_tranchSNP -truth-sensitivity-filter-level 99.8 --create-output-variant-index true -mode SNP
# Apply recalibration to Indels/Apply Indel VQSR on genotypes file (INDEL sensitivity 99.0)/Apply Indel VQSR on genotypes file
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=20" ApplyVQSR -V $recal_snpsindel -R $ref -L $hg38exon_interval -O $vcfFile_recal_v2 --recal-file $output_recalINDEL --tranches-file $output_tranchINDEL -truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL

exit 0

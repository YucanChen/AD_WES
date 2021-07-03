#!/bin/bash

GATK=~/gatk-4.1.0.0/gatk
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
hg38exon_interval=~/gatk-4.1.0.0/bundle_hg38/targets.preprocessed.interval.list

vcfHapmap=~/gatk-4.1.0.0/bundle_hg38/hapmap_3.3.hg38.vcf.gz
vcfOmni=~/gatk-4.1.0.0/bundle_hg38/1000G_omni2.5.hg38.vcf.gz
vcfGlk=~/gatk-4.1.0.0/bundle_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
vcfDbsnp=~/gatk-4.1.0.0/bundle_hg38/dbsnp_146.hg38.vcf.gz
vcfMills=~/gatk-4.1.0.0/bundle_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

vcfFile=/WES_output_combine/combine_genotype_vcf/WES_combine.vcf
vcfFile_sites_only=/WES_output_combine/combine_genotype_vcf/WES_combine_sites_only.vcf.gz

output_recalSNP=/WES_output_combine/VQSR/WES_recalibrate_SNP.recal
output_tranchSNP=/WES_output_combine/VQSR/WES_recalibrate_SNP.tranches
output_rscriptSNP=/WES_output_combine/VQSR/output.SNP.recalibration_plots.rscript
recal_snpsindel=/WES_output_combine/VQSR/WES_recalibrated_snps_raw_indels.vcf

output_recalINDEL=/WES_output_combine/VQSR/WES_recalibrate_INDEL.recal
output_tranchINDEL=/WES_output_combine/VQSR/WES_recalibrate_INDEL.tranches
output_rscriptINDEL=/WES_output_combine/VQSR/output.INDEL.recalibration_plots.rscript
vcfFile_recal=/WES_output_combine/VQSR/WES_combine.recal.vcf

# Strip to sites only
$GATK MakeSitesOnlyVcf --INPUT $vcfFile --OUTPUT $vcfFile_sites_only

# VARIANT QUALITY SCORE RECALIBRATION - SNPs
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=40" VariantRecalibrator -V $vcfFile_sites_only -R $ref -O $output_recalSNP -mode SNP --tranches-file $output_tranchSNP -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an FS -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -an MQ -an SOR --max-gaussians 6 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $vcfHapmap -resource:omni,known=false,training=true,truth=true,prior=12.0 $vcfOmni -resource:1000G,known=false,training=true,truth=false,prior=10.0 $vcfGlk -resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $vcfDbsnp --rscript-file $output_rscriptSNP

# Run Variant Recalibrator – Indels
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=40" VariantRecalibrator -V $vcfFile_sites_only -R $ref -O $output_recalINDEL -mode INDEL --tranches-file $output_tranchINDEL -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an InbreedingCoeff -an MQRankSum -an QD -an SOR --max-gaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 $vcfMills -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $vcfDbsnp --rscript-file $output_rscriptINDEL
wait

# Apply recalibration to SNPs (SNP sensitivity 99.6)
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=20" ApplyVQSR -V $vcfFile -R $ref -L $hg38exon_interval -O $recal_snpsindel  --recal-file $output_recalSNP --tranches-file $output_tranchSNP -truth-sensitivity-filter-level 99.6 --create-output-variant-index true -mode SNP

# Apply recalibration to Indels/Apply Indel VQSR on genotypes file (INDEL sensitivity 95.0)
$GATK --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=20" ApplyVQSR -V $recal_snpsindel -R $ref -L $hg38exon_interval -O $vcfFile_recal --recal-file $output_recalINDEL --tranches-file $output_tranchINDEL -truth-sensitivity-filter-level 95.0 --create-output-variant-index true -mode INDEL

exit 0

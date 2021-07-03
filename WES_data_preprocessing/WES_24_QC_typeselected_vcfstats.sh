#!/bin/bash

snp=/WES_output_combine/QC_passed_vcf/WES_type_selected.snp.vcf
indel=/WES_output_combine/QC_passed_vcf/WES_type_selected.indel.vcf
outstats_snp=/WES_output_combine/QC_passed_vcf/result_vcfstats_typeselected/WES_highquality.snpstats
outstats_indel=/WES_output_combine/QC_passed_vcf/result_vcfstats_typeselected/WES_highquality.indelstats
rtg vcfstats $snp > $outstats_snp &
rtg vcfstats $indel > $outstats_indel &
wait

exit 0

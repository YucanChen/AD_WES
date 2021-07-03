#!/bin/bash

cd /afterrecall_QC/vcftools_relatedness/
vcftools --relatedness2 --vcf /WES_output_combine/VQSR/WES_combine.recal_v2.vcf --out WES_combine.recal_v2

exit 0

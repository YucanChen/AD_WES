#!/bin/bash

input_snp=/WES_output_combine/QC_passed_vcf/WES_type_selected.snp.vcf
input_indel=/WES_output_combine/QC_passed_vcf/WES_type_selected.indel.vcf
output_indel=/WES_output_combine/ANNOVAR_annotation/WES_highquality.indel
output_snp=/WES_output_combine/ANNOVAR_annotation/WES_highquality.snp

cd ~/biosoft/annovar/
table_annovar.pl $input_indel humandb/ -buildver hg38 -out $output_indel -remove -protocol refGene,cytoBand,dbscsnv11,avsnp147,exac03,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,gnomad211_exome,X1000G_AFR,X1000G_ALL,X1000G_AMR,X1000G_EAS,X1000G_EUR,X1000G_SAS,nci60,cosmic70,clinvar_20200316,dbnsfp33a -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish &
table_annovar.pl $input_snp humandb/ -buildver hg38 -out $output_snp -remove -protocol refGene,cytoBand,dbscsnv11,avsnp147,exac03,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,gnomad211_exome,X1000G_AFR,X1000G_ALL,X1000G_AMR,X1000G_EAS,X1000G_EUR,X1000G_SAS,nci60,cosmic70,clinvar_20200316,dbnsfp33a -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish &
wait

exit 0

#!/bin/bash

snpsubset_AD=/WES_output_combine/QC_passed_vcf/WES_type_selected_AD.snp.vcf
indelsubset_AD=/WES_output_combine/QC_passed_vcf/WES_type_selected_AD.indel.vcf
snpsubset_ADearly=/WES_output_combine/QC_passed_vcf/WES_type_selected_ADearly.snp.vcf
indelsubset_ADearly=/WES_output_combine/QC_passed_vcf/WES_type_selected_ADearly.indel.vcf
output_indel_AD=/WES_output_combine/ANNOVAR_annotation/WES_highquality_AD.indel
output_snp_AD=/WES_output_combine/ANNOVAR_annotation/WES_highquality_AD.snp
output_indel_ADearly=/WES_output_combine/ANNOVAR_annotation/WES_highquality_ADearly.indel
output_snp_ADearly=/WES_output_combine/ANNOVAR_annotation/WES_highquality_ADearly.snp

cd ~/biosoft/annovar/
table_annovar.pl $indelsubset_AD humandb/ -buildver hg38 -out $output_indel_AD -remove -protocol refGene,cytoBand,dbscsnv11,avsnp147,exac03,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,gnomad211_exome,X1000G_AFR,X1000G_ALL,X1000G_AMR,X1000G_EAS,X1000G_EUR,X1000G_SAS,nci60,cosmic70,clinvar_20200316,dbnsfp33a -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish &
table_annovar.pl $snpsubset_AD humandb/ -buildver hg38 -out $output_snp_AD -remove -protocol refGene,cytoBand,dbscsnv11,avsnp147,exac03,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,gnomad211_exome,X1000G_AFR,X1000G_ALL,X1000G_AMR,X1000G_EAS,X1000G_EUR,X1000G_SAS,nci60,cosmic70,clinvar_20200316,dbnsfp33a -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish &

table_annovar.pl $indelsubset_ADearly humandb/ -buildver hg38 -out $output_indel_ADearly -remove -protocol refGene,cytoBand,dbscsnv11,avsnp147,exac03,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,gnomad211_exome,X1000G_AFR,X1000G_ALL,X1000G_AMR,X1000G_EAS,X1000G_EUR,X1000G_SAS,nci60,cosmic70,clinvar_20200316,dbnsfp33a -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish &
table_annovar.pl $snpsubset_ADearly humandb/ -buildver hg38 -out $output_snp_ADearly -remove -protocol refGene,cytoBand,dbscsnv11,avsnp147,exac03,esp6500siv2_ea,esp6500siv2_aa,esp6500siv2_all,gnomad211_exome,X1000G_AFR,X1000G_ALL,X1000G_AMR,X1000G_EAS,X1000G_EUR,X1000G_SAS,nci60,cosmic70,clinvar_20200316,dbnsfp33a -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish &
wait

exit 0

#!/bin/bash
cd /gvcf_filtered/gvcf_depth/
for i in /gvcf/*_dedup_recal.g.vcf
do
{
j=${i##*/gvcf/}
vcftools --vcf $i --depth --out /gvcf_filtered/gvcf_depth/$j
average_depth=$(awk 'NR==2{print $3}' $j.idepth)
depth_cutoff=$(echo "$average_depth + 4 * sqrt($average_depth)"|bc)
qual_cutoff=$(echo "$depth_cutoff*2.0"|bc)
tmp1=/gvcf_filtered/gvcf_tmp/$j.tmp1
tmp2=/gvcf_filtered/gvcf_tmp/$j.tmp2
tmp3=/gvcf_filtered/gvcf_tmp/$j.tmp3
tmp4=/gvcf_filtered/gvcf_tmp/$j.tmp4
final_vcf=/gvcf_filtered/${j%_dedup_recal.g.vcf}_dedup_recal_highDPrm.g.vcf
vcftools --vcf $i --min-meanDP $depth_cutoff --stdout --recode > $tmp1 && bcftools filter -e 'GT="RR"' $tmp1 > $tmp2 && vcftools --vcf $tmp2 --minQ $qual_cutoff --stdout --recode > $tmp3 && vcftools --vcf $tmp2 --exclude-positions $tmp3 --stdout --recode > $tmp4 && vcftools --vcf $i --exclude-positions $tmp4 --stdout --recode --recode-INFO-all > $final_vcf && rm $tmp1 && rm $tmp2 && rm $tmp3
}&
done
wait

exit 0

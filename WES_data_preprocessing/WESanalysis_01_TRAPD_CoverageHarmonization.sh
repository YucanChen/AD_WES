#!/bin/bash

# section 1
ref=~/gatk-4.1.0.0/bundle_hg38/Homo_sapiens_assembly38.fasta
output=~/burden_test/case/cases.counts_samtools.txt
bam_list=~/burden_test/case/case_bams.list
hg38exon_bed=~/bed/hg38.exon.bed

samtools depth -b $hg38exon_bed -f $bam_list -q 0 -Q 20 --reference $ref > $output

# section 2
case_count=~/burden_test/case/cases.counts_samtools.txt
case_bed=~/burden_test/case/cases.dp10.bed
cat $case_count | awk '{count=0} {for(i=3; i<1000; i++) if($i>10) count++} {if(count>0.9) print $1"\t"$2}' | awk '{print $1"\t"($2-1)"\t"$2}' | bedtools merge -i stdin > $case_bed

# section 3
control_bed=~/burden_test/control_data/gnomad.dp10.map.bed
case_bed=~/burden_test/case/cases.dp10.bed
combine_bed=~/burden_test/case/combined.dp10.bed
bedtools intersect -a $control_bed -b $case_bed | sort -k1.4,1n -k2,2n | bedtools merge -i stdin > $combine_bed

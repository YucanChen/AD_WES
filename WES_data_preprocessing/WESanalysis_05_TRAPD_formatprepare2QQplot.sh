#### make_snp_file
bedfile=~/burden_test/case/combined.dp10.bed

vcffile_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.vcf.gz
outfile_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.filter.snpfile.txt

vcffile_case_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.vcf.gz
outfile_case_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.filter.indelfile.txt

vcffile_control=~/burden_test/control_data/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.vcf.gz
outfile_control=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.filter.snpfile.txt
vcffile_control_indel=~/burden_test/control_data/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTV_indel2.vcf.gz
outfile_control_indel=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTV.filter.indelfile.txt

vcffile_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.vcf.bgz
outfile_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.filter.snpfile.txt

vcffile_case_indel_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.vcf.bgz
outfile_case_indel_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.filter.indelfile.txt

python ~/TRAPD/code/make_snp_file.py --vcffile $vcffile_case --outfile $outfile_case --genecolname "Gene.refGene" --bedfile $bedfile --snponly &
python ~/TRAPD/code/make_snp_file.py --vcffile $vcffile_case_indel --outfile $outfile_case_indel --genecolname "Gene.refGene" --bedfile $bedfile --indelonly &
python ~/TRAPD/code/make_snp_file.py --vcffile $vcffile_control --outfile $outfile_control --vep --genecolname "SYMBOL" --pass --snponly --bedfile $bedfile --includeinfo "AF_popmax[<]0.001" --includeinfo "non_topmed_AF_popmax[<]0.001" --includeinfo "non_neuro_AF_popmax[<]0.001" --includeinfo "non_cancer_AF_popmax[<]0.001" --includeinfo "controls_AF_popmax[<]0.001" &
python ~/TRAPD/code/make_snp_file.py --vcffile $vcffile_control_indel --outfile $outfile_control_indel --vep --genecolname "SYMBOL" --pass --indelonly --bedfile $bedfile --includeinfo "AF_popmax[<]0.001" --includeinfo "non_topmed_AF_popmax[<]0.001" --includeinfo "non_neuro_AF_popmax[<]0.001" --includeinfo "non_cancer_AF_popmax[<]0.001" --includeinfo "controls_AF_popmax[<]0.001" &

python ~/TRAPD/code/make_snp_file.py --vcffile $vcffile_case_inner --outfile $outfile_case_inner --genecolname "Gene.refGene" --bedfile $bedfile --snponly &
python ~/TRAPD/code/make_snp_file.py --vcffile $vcffile_case_indel_inner --outfile $outfile_case_indel_inner --genecolname "Gene.refGene" --bedfile $bedfile --indelonly &
wait

#### 1b) Merging SNP files This is an optional step if you need to merge two SNP files (for example, if you performed step 1a separately for SNPs and indels). It can also be used if you perfomed Step 1a separately for each chromosome. It has two required options:
snpfile_case=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.filter.snpfile.txt
indelfile_case=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY.filter.indelfile.txt
snpfile_control=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.filter.snpfile.txt
indelfile_control=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTV.filter.indelfile.txt
outfile_case=~/burden_test/case/WES_highquality_hg38_sorted.all_deXY_PAVcal.filter.outfile.txt
outfile_control=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTVcal.filter.outfile.txt

snpfile_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.filter.snpfile.txt
indelfile_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.filter.indelfile.txt
outfile_case_inner=~/burden_test/case/WES_highquality_hg38_sorted.all_deXY_PAVcal_inmaf.filter.outfile.txt

python ~/TRAPD/code/merge_snp_file.py --snpfiles $snpfile_case,$indelfile_case --outfile $outfile_case &
python ~/TRAPD/code/merge_snp_file.py --snpfiles $snpfile_control,$indelfile_control --outfile $outfile_control &
python ~/TRAPD/code/merge_snp_file.py --snpfiles $snpfile_case_inner,$indelfile_case_inner --outfile $outfile_case_inner &
wait

#### 2a) Counting carriers in case cohort This script will tabulate the number of cases carrying qualifying variants in each gene as defined by a SNP file.
#The command takes in a vcf file containing case sample genotypes and a SNP file listing the qualifying variants for each gene. The general command is:
bedfile=~/burden_test/case/combined.dp10.bed
case_vcf=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.vcf.gz
case_snpfile=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.filter.snpfile.txt
casecounts=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.filter.casecounts.txt

case_vcf_merge=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY.merge.vcf.gz
case_snpfile_merge=~/burden_test/case/WES_highquality_hg38_sorted.all_deXY_PAVcal.filter.outfile.txt
casecounts_merge=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY.filter.casecounts.merge.txt

case_vcf_merge_inner=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY_inmaf.merge.vcf.gz
case_snpfile_merge_inner=~/burden_test/case/WES_highquality_hg38_sorted.all_deXY_PAVcal_inmaf.filter.outfile.txt
casecounts_merge_inner=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY_inmaf.filter.casecounts.merge.txt

python ~/TRAPD/code/count_cases.py -v $case_vcf -s $case_snpfile -o $casecounts --pass --bedfile $bedfile &
python ~/TRAPD/code/count_cases.py -v $case_vcf_merge -s $case_snpfile_merge -o $casecounts_merge --pass --bedfile $bedfile &
python ~/TRAPD/code/count_cases.py -v $case_vcf_merge_inner -s $case_snpfile_merge_inner -o $casecounts_merge_inner --pass --bedfile $bedfile &
wait

#### 2b) Counting carriers in public control cohorts
#This script will tabulate the approximate number of controls carrying qualifying variants in each gene as defined by a SNP file. Currently, this script has been configured to run using ExAC (http://exac.broadinstitute.org/downloads) or gnomAD (http://gnomad.broadinstitute.org/) data, but any sites-only vcf can be used (as long as it contains AC and AN in the INFO field).
#The general command is:
bedfile=~/burden_test/case/combined.dp10.bed
control_vcf=~/burden_test/control_data/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.vcf.gz
control_snpfile=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.filter.snpfile.txt
controlcounts=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.filter.casecounts.txt

control_vcf_merge=~/burden_test/control_data/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTVcal.merge.vcf.gz
control_snpfile_merge=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTVcal.filter.outfile.txt
controlcounts_merge=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTVcal.filter.casecounts.merge.txt

python ~/TRAPD/code/count_controls.py -v $control_vcf -s $control_snpfile -o $controlcounts --pass --bedfile $bedfile --database "gnomad"
python ~/TRAPD/code/count_controls.py -v $control_vcf_merge -s $control_snpfile_merge -o $controlcounts_merge --pass --bedfile $bedfile --database "gnomad"

#### 3) Run burden testing This script will run the actual burden testing. It performs a one-sided Fisher's exact test to determine if there is a greater burden of qualifying variants in cases as compared to controls for each gene. It will perform this burden testing under a dominant and a recessive model.
#It requires R; the script was tested using R v3.1, but any version of R should work. The script should be run as:
casecounts=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY.filter.casecounts.txt
controlcounts=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_D_PTV.filter.casecounts.txt
outfile=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_D_PTV.txt

casecounts_merge=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY.filter.casecounts.merge.txt
controlcounts_merge=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTVcal.filter.casecounts.merge.txt
outfile_merge=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_PTVcal.merge.txt

casecounts_merge_inner=~/burden_test/case/WES_highquality_hg38_sorted.PAV_D_deXY_inmaf.filter.casecounts.merge.txt
controlcounts_merge=~/burden_test/case/gnomad.exomes.r2.1.1.sites.liftover_grch38_deXYclean_PTVcal.filter.casecounts.merge.txt
outfile_merge_inner=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_PTVcal_inmaf.merge.txt

Rscript=~/R-2.15.3/bin/Rscript
burden_R=~/TRAPD/code/burden.R

$Rscript $burden_R --casefile $casecounts --casesize 60 --controlfile $controlcounts --controlsize 125748 --outfile $outfile
$Rscript $burden_R --casefile $casecounts_merge --casesize 60 --controlfile $controlcounts_merge --controlsize 125748 --outfile $outfile_merge
$Rscript $burden_R --casefile $casecounts_merge_inner --casesize 60 --controlfile $controlcounts_merge --controlsize 125748 --outfile $outfile_merge_inner

#### Generate QQ Plot The last step is the generate the QQ plot and calculate the lambda_delta95 as described in our paper
burden_outfile=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_D_PTV.txt
out_png=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_D_PTV.png
burden_outfile_merge=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_PTVcal.merge.txt
out_png_merge=~/burden_test/burden_test_results/burden.out_grch38_deXYclean_D_PTV.merge.png
Rscript=~/R-2.15.3/bin/Rscript

$Rscript ~/TRAPD/code/QQ.R --pvalfile $burden_outfile --plotfile $out_png --removenull &
$Rscript ~/TRAPD/code/QQ.R --pvalfile $burden_outfile_merge --plotfile $out_png_merge --removenull &
wait

#### vcf2bed_snpindel
vcf_snp=~/burden_test/case/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.vcf.bgz
vcf_indel=~/burden_test/case/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.vcf.bgz

tmp_snp=~/burden_test/burden_test_results/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.vcf
tmp_indel=~/burden_test/burden_test_results/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.vcf

bed_snp=~/burden_test/burden_test_results/WES_highquality_hg38_sorted.snp.PAV_D_deXY_inmaf.bed
bed_indel=~/burden_test/burden_test_results/WES_highquality_hg38_sorted.indel.PAV_deXY_inmaf.bed

zcat $vcf_snp > $tmp_snp && convert2bed --input=vcf --output=bed < $tmp_snp > $bed_snp &
zcat $vcf_indel > $tmp_indel && convert2bed --input=vcf --output=bed < $tmp_indel > $bed_indel &
wait

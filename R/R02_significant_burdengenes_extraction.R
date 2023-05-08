library(openxlsx)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))

setwd("/AD_exome/Downstream_analyses/AD_downstream_scripts_20220627/")
outPath<-'/AD_exome/Downstream_analyses/AD_downstream_scripts_20220627/burden_test'

## import WES burden genes:
burdengenes_dementia=read.table("/AD_exome/CYC_recall_batch2/WES_analysis_3rdadjust/burden_test/burden_test_results_dementia/burden.out_grch38_deXYclean_PTVcal_inmaf.merge.txt",header = T)
burdengenes_AD=read.table("/AD_exome/CYC_recall_batch2/WES_analysis_3rdadjust/burden_test/burden_test_results_AD/burden.out_grch38_deXYclean_PTVcal_inmaf.merge.txt",header = T)

## genename format transformation:
burdengenes_dementia.names <- gsub("\\\\x3b.*","",burdengenes_dementia$GENE)
burdengenes_dementia.names[1:100]
burdengenes_dementia.subnames <- gsub(".*\\\\x3b","",burdengenes_dementia$GENE)
burdengenes_dementia.subnames[1:100]
burdengenes_dementia$GENE <- burdengenes_dementia.names
burdengenes_dementia$GENE_sub <- burdengenes_dementia.subnames

burdengenes_AD.names <- gsub("\\\\x3b.*","",burdengenes_AD$GENE)
burdengenes_AD.names[1:100]
burdengenes_AD.subnames <- gsub(".*\\\\x3b","",burdengenes_AD$GENE)
burdengenes_AD.subnames[1:100]
burdengenes_AD$GENE <- burdengenes_AD.names
burdengenes_AD$GENE_sub <- burdengenes_AD.subnames


## significant burden genes selection:
burdengenes_dementia$P_DOM_adjustFDR<-p.adjust(burdengenes_dementia$P_DOM,"fdr")
burdengenes_AD$P_DOM_adjustFDR<-p.adjust(burdengenes_AD$P_DOM,"fdr")

burdengenes_dementia$P_DOM_adjustbonf<-p.adjust(burdengenes_dementia$P_DOM,"bonferroni")
burdengenes_AD$P_DOM_adjustbonf<-p.adjust(burdengenes_AD$P_DOM,"bonferroni")

burdengenes_dementia$exist_mut <- burdengenes_dementia$CASE_TOTAL_AC!=0
burdengenes_AD$exist_mut<-burdengenes_AD$CASE_TOTAL_AC!=0


## output all results:
xl_list <- list()
xl_list$dementia <- burdengenes_dementia %>% dplyr::arrange(P_DOM_adjustFDR)
xl_list$AD <- burdengenes_AD %>% dplyr::arrange(P_DOM_adjustFDR)
write.xlsx(xl_list,paste(outPath,'burdengenes_result_extraction.xlsx',sep = '/'))

burdengenes_AD_list <- burdengenes_AD[which(burdengenes_AD$P_DOM_adjustFDR < 0.05),c("GENE","GENE_sub")]
burdengenes_dementia_list <- burdengenes_dementia[which(burdengenes_dementia$P_DOM_adjustFDR < 0.05),c("GENE","GENE_sub")]

write.csv(burdengenes_dementia_list,paste(outPath,"list_burdengenes_dementia.csv",sep = '/'),row.names = FALSE,quote = F)
write.csv(burdengenes_AD_list,paste(outPath,"list_burdengenes_AD.csv",sep = '/'),row.names = FALSE,quote = F)

## significant gene list storage:
save(burdengenes_AD_list,burdengenes_dementia_list,file = paste(str_replace(outPath,'/burden_test',''),"part1_0_significant_burdengenes_FDR05.RData",sep = '/'))

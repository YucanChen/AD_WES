library(openxlsx)
setwd("/AD_exome/Downstream_analyses/AD_downstream_scripts/")

###import WES burden genes:
burdengenes_dementia=read.xlsx("/AD_exome/Downstream_analyses/part1_burdentest/burden_test_batch2_allresults/burden.out_grch38_deXYclean_PTVcal_inmaf_merge_all.xlsx",sheet = "dementia")
burdengenes_AD=read.xlsx("/AD_exome/Downstream_analyses/part1_burdentest/burden_test_batch2_allresults/burden.out_grch38_deXYclean_PTVcal_inmaf_merge_all.xlsx",sheet = "AD")

###genename format transformation:
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

###significant burden genes selection:
burdengenes_dementia$P_DOM_adjustFDR<-p.adjust(burdengenes_dementia$P_DOM,"fdr")
burdengenes_AD$P_DOM_adjustFDR<-p.adjust(burdengenes_AD$P_DOM,"fdr")

burdengenes_dementia$P_DOM_adjustbonf<-p.adjust(burdengenes_dementia$P_DOM,"bonferroni")
burdengenes_AD$P_DOM_adjustbonf<-p.adjust(burdengenes_AD$P_DOM,"bonferroni")

burdengenes_dementia$exist_mut <- burdengenes_dementia$CASE_TOTAL_AC!=0
burdengenes_AD$exist_mut<-burdengenes_AD$CASE_TOTAL_AC!=0


###output all results:

xl_list <- list()
xl_list$dementia <- burdengenes_dementia
xl_list$AD <- burdengenes_AD

write.xlsx(xl_list,"F:/Lin_Lab/AD_exome/Downstream_analyses/part1_burdentest/burden_test_batch2_allresults/burdengenes_result_extraction.xlsx")
write.xlsx(xl_list,"burdengenes_result_extraction.xlsx")

burdengenes_AD_list <- burdengenes_AD[which(burdengenes_AD$P_DOM_adjustFDR < 0.05),c("GENE","GENE_sub")]
burdengenes_dementia_list <- burdengenes_dementia[which(burdengenes_dementia$P_DOM_adjustFDR < 0.05),c("GENE","GENE_sub")]

write.csv(burdengenes_dementia_list,"list_burdengenes_dementia.csv",row.names = FALSE)
write.csv(burdengenes_AD_list,"list_burdengenes_AD.csv",row.names = FALSE)

###significant gene list storage:
save(burdengenes_AD_list,burdengenes_ADearly_list,burdengenes_dementia_list,file = "part1_0_significant_burdengenes_FDR05.RData")

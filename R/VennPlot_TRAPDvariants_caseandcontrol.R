library(openxlsx)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)

# Venn plot
## dementia
Venn.case_dementia <- read.xlsx("case_allgroup_outfile_snps.xlsx",sheet = "case_dementia")[,"GENE"]
Venn.control_dementia <- read.table("control_dementia_outfile_snps.txt",sep = '\t',header = T)[,"GENE"]
list_combined_dementia <- list("case variants"=Venn.case_dementia,"control variants"=Venn.control_dementia)
venn.diagram(list_combined_dementia,height = 3000, width = 3000,cat.dist= 0.25,cat.pos=c(0,0),
             resolution = 300, imagetype = "png", alpha=c(0.45,0.45),cex=2.5,
             fill=c("blue","yellow"), cat.fontface=4,fontfamily=3,cat.cex=3,sub.cex=3,
             filename = "../AD_paperwrittingplots/VennDiagram_TRAPDvariants_casecontrol_dementia.png")
Venndata_dementia <- get.venn.partitions(list_combined_dementia)
for (i in 1:nrow(Venndata_dementia)){
  Venndata_dementia[i,'values'] <- paste(Venndata_dementia[[i,'..values..']], collapse = ', ')
}
write.table(Venndata_dementia[-c(4)], 'part5_3_Venndata_dementia.txt', row.names = FALSE, sep = '\t', quote = FALSE)
Venndata_dementia.caseonly <- Venndata_dementia[3,"..values.."]

## AD
Venn.case_AD <- read.xlsx("case_allgroup_outfile_snps.xlsx",sheet = "case_AD")[,"GENE"]
Venn.control_AD <- read.table("control_AD_outfile_snps.txt",sep = '\t',header = T)[,"GENE"]
list_combined_AD <- list("case variants"=Venn.case_AD,"control variants"=Venn.control_AD)
venn.diagram(list_combined_AD,height = 3000, width = 3000,cat.dist= 0.25,cat.pos=c(0,0),
             resolution = 300, imagetype = "png", alpha=c(0.45,0.45),cex=2.5,
             fill=c("blue","yellow"), cat.fontface=4,fontfamily=3,cat.cex=3,sub.cex=3,
             filename = "../AD_paperwrittingplots/VennDiagram_TRAPDvariants_casecontrol_AD.png")
Venndata_AD <- get.venn.partitions(list_combined_AD)
for (i in 1:nrow(Venndata_AD)){
  Venndata_AD[i,'values'] <- paste(Venndata_AD[[i,'..values..']], collapse = ', ')
}
write.table(Venndata_AD[-c(4)], 'part5_3_Venndata_AD.txt', row.names = FALSE, sep = '\t', quote = FALSE)
Venndata_AD.caseonly <- Venndata_AD[3,"..values.."]

## early AD
Venn.case_earlyAD <- read.xlsx("case_allgroup_outfile_snps.xlsx",sheet = "case_earlyAD")[,"GENE"]
Venn.control_earlyAD <- read.table("control_earlyAD_outfile_snps.txt",sep = '\t',header = T)[,"GENE"]
list_combined_earlyAD <- list("case variants"=Venn.case_earlyAD,"control variants"=Venn.control_earlyAD)
venn.diagram(list_combined_earlyAD,height = 3000, width = 3000,cat.dist= 0.25,cat.pos=c(0,0),
             resolution = 300, imagetype = "png", alpha=c(0.45,0.45),cex=2.5,
             fill=c("blue","yellow"), cat.fontface=4,fontfamily=3,cat.cex=3,sub.cex=3,
             filename = "../AD_paperwrittingplots/VennDiagram_TRAPDvariants_casecontrol_earlyAD.png")
Venndata_earlyAD <- get.venn.partitions(list_combined_earlyAD)
for (i in 1:nrow(Venndata_earlyAD)){
  Venndata_earlyAD[i,'values'] <- paste(Venndata_earlyAD[[i,'..values..']], collapse = ', ')
}
write.table(Venndata_earlyAD[-c(4)], 'part5_3_Venndata_earlyAD.txt', row.names = FALSE, sep = '\t', quote = FALSE)
Venndata_earlyAD.caseonly <- Venndata_earlyAD[3,"..values.."]

# Bar plot: mutation frequency
freq_case <- as.data.frame(table(Venn.case_dementia));colnames(freq_case) <- c("GENE","case_freq")
freq_control <- as.data.frame(table(Venn.control_dementia));colnames(freq_control) <- c("GENE","case_control")
mutationfreq_dementia <- merge(freq_case,freq_control,by="GENE");mutationfreq_dementia <- mutationfreq_dementia[order(mutationfreq_dementia$case_freq,decreasing = T),]
casefreq_dementia <- as.data.frame(table(mutationfreq_dementia$case_freq));colnames(casefreq_dementia) <- c("mutation_freq","gene_number")
casefreq_dementia$group <- rep("dementia",length(casefreq_dementia$mutation_freq))

freq_case <- as.data.frame(table(Venn.case_AD));colnames(freq_case) <- c("GENE","case_freq")
freq_control <- as.data.frame(table(Venn.control_AD));colnames(freq_control) <- c("GENE","case_control")
mutationfreq_AD <- merge(freq_case,freq_control,by="GENE");mutationfreq_AD <- mutationfreq_AD[order(mutationfreq_AD$case_freq,decreasing = T),]
casefreq_AD <- as.data.frame(table(mutationfreq_AD$case_freq));colnames(casefreq_AD) <- c("mutation_freq","gene_number")
casefreq_AD$group <- rep("AD",length(casefreq_AD$mutation_freq))

freq_case <- as.data.frame(table(Venn.case_earlyAD));colnames(freq_case) <- c("GENE","case_freq")
freq_control <- as.data.frame(table(Venn.control_earlyAD));colnames(freq_control) <- c("GENE","case_control")
mutationfreq_earlyAD <- merge(freq_case,freq_control,by="GENE");mutationfreq_earlyAD <- mutationfreq_earlyAD[order(mutationfreq_earlyAD$case_freq,decreasing = T),]
casefreq_earlyAD <- as.data.frame(table(mutationfreq_earlyAD$case_freq));colnames(casefreq_earlyAD) <- c("mutation_freq","gene_number")
casefreq_earlyAD$group <- rep("earlyAD",length(casefreq_earlyAD$mutation_freq))

rm(Venn.control_AD,Venn.control_dementia,Venn.control_earlyAD)
combine_casefreq <- rbind(casefreq_earlyAD,casefreq_AD,casefreq_dementia)

ggplot(data=combine_casefreq)+
  geom_bar(aes(x=as.numeric(as.character(mutation_freq)),y=gene_number,fill=group),stat='identity',position = "dodge",width=0.6) +
  scale_fill_manual(values = c("#990000","#0066CC","#660099"))+
  coord_flip() +
  xlab("mutation frequency") +
  ylab("gene number") +
  theme(
    axis.text.x=element_text(color="black",size=rel(1.1)),
    axis.text.y=element_text(color="black", size=rel(1.1)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    legend.text=element_text(color="black",size=rel(0.9)),
    legend.title = element_text(color="black",size=rel(1.0))
  )

# Functional enrichment analysis
## clusterProfiler enrichment analysis
ID_conversion_dementia <- read.csv("E:/Lin_Lab/AD_exome/Downstream_analyses/AD_downstream_scripts/burdenresult_dementia_annotatedforLR.csv")[,c("entrezgene_id","GENE")]
ID_conversion_AD <- read.csv("E:/Lin_Lab/AD_exome/Downstream_analyses/AD_downstream_scripts/burdenresult_AD_annotatedforLR.csv")[,c("entrezgene_id","GENE")]
ID_conversion_earlyAD <- read.csv("E:/Lin_Lab/AD_exome/Downstream_analyses/AD_downstream_scripts/burdenresult_earlyAD_annotatedforLR.csv")[,c("entrezgene_id","GENE")]

Venndata_dementia.caseonly_entrezid <- unlist(lapply(Venndata_dementia.caseonly, function(x){ID_conversion_dementia[which(ID_conversion_dementia$GENE %in% x),"entrezgene_id"]}))
Venndata_AD.caseonly_entrezid <- unlist(lapply(Venndata_AD.caseonly, function(x){ID_conversion_AD[which(ID_conversion_AD$GENE %in% x),"entrezgene_id"]}))
Venndata_earlyAD.caseonly_entrezid <- unlist(lapply(Venndata_earlyAD.caseonly, function(x){ID_conversion_earlyAD[which(ID_conversion_earlyAD$GENE %in% x),"entrezgene_id"]}))

enrichBP_caseonly_dementia <- enrichGO(Venndata_dementia.caseonly_entrezid,OrgDb="org.Hs.eg.db",keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.1,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 10,maxGSSize = 1000,readable = FALSE,pool = FALSE)
enrichBP_caseonly_dementia <- as.data.frame(enrichBP_caseonly_dementia)
enrichBP_caseonly_AD <- enrichGO(Venndata_AD.caseonly_entrezid,OrgDb="org.Hs.eg.db",keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.1,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 10,maxGSSize = 1000,readable = FALSE,pool = FALSE)
enrichBP_caseonly_AD <- as.data.frame(enrichBP_caseonly_AD)
enrichBP_caseonly_earlyAD <- enrichGO(Venndata_earlyAD.caseonly_entrezid,OrgDb="org.Hs.eg.db",keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.1,pAdjustMethod = "BH",qvalueCutoff = 0.1,minGSSize = 10,maxGSSize = 1000,readable = FALSE,pool = FALSE)
enrichBP_caseonly_earlyAD <- as.data.frame(enrichBP_caseonly_earlyAD)

sheets = list("enrichBP_dementia"=enrichBP_caseonly_dementia,"enrichBP_AD"=enrichBP_caseonly_AD,"enrichBP_earlyAD"=enrichBP_caseonly_earlyAD)
write.xlsx(sheets,"part5_3_enrichBP_caseonly_allgroup_result_BH1.xlsx")

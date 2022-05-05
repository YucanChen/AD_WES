library(openxlsx)
library(ggplot2)
library(VennDiagram)
library(devtools)
library(Vennerable)


# Venn

## dementia

Venn.case_dementia <- read.xlsx("./data/burdenresult_dementia/case_allgroup_outfile_snps.xlsx",sheet = "case_dementia")[,"GENE"]
Venn.control_dementia <- read.table("./data/burdenresult_dementia/control_dementia_outfile_snps.txt",sep = '\t',header = T)[,"GENE"]
list_combined_dementia <- list("case variants"=Venn.case_dementia,"control variants"=Venn.control_dementia)
Venndata_dementia <- get.venn.partitions(list_combined_dementia)
for (i in 1:nrow(Venndata_dementia)){
  Venndata_dementia[i,'values'] <- paste(Venndata_dementia[[i,'..values..']], collapse = ', ')
}
Venndata_dementia.caseonly <- Venndata_dementia[3,"..values.."]

## AD

Venn.case_AD <- read.xlsx("./data/burdenresult_AD/case_allgroup_outfile_snps.xlsx",sheet = "case_AD")[,"GENE"]
Venn.control_AD <- read.table("./data/burdenresult_AD/control_AD_outfile_snps.txt",sep = '\t',header = T)[,"GENE"]
list_combined_AD <- list("case variants"=Venn.case_AD,"control variants"=Venn.control_AD)
Venndata_AD <- get.venn.partitions(list_combined_AD)
for (i in 1:nrow(Venndata_AD)){
  Venndata_AD[i,'values'] <- paste(Venndata_AD[[i,'..values..']], collapse = ', ')
}
Venndata_AD.caseonly <- Venndata_AD[3,"..values.."]

## early AD

Venn.case_earlyAD <- read.xlsx("./data/burdenresult_earlyAD/case_allgroup_outfile_snps.xlsx",sheet = "case_earlyAD")[,"GENE"]
Venn.control_earlyAD <- read.table("./data/burdenresult_earlyAD/control_earlyAD_outfile_snps.txt",sep = '\t',header = T)[,"GENE"]
list_combined_earlyAD <- list("case variants"=Venn.case_earlyAD,"control variants"=Venn.control_earlyAD)
Venndata_earlyAD <- get.venn.partitions(list_combined_earlyAD)
for (i in 1:nrow(Venndata_earlyAD)){
  Venndata_earlyAD[i,'values'] <- paste(Venndata_earlyAD[[i,'..values..']], collapse = ', ')
}
Venndata_earlyAD.caseonly <- Venndata_earlyAD[3,"..values.."]

# Entrez id mapping and output results

ID_mapping<-read.csv("./metadata/burdenresult_dementia_annotatedforLR.csv")[,c(1,2,5)]

ENTREZID_case_dementia<-unique(unlist(lapply(Venn.case_dementia,function(x){ID_mapping[which(ID_mapping$GENE==x),"entrez_id"]})))
ENTREZID_case_AD<-unique(unlist(lapply(Venn.case_AD,function(x){ID_mapping[which(ID_mapping$GENE==x),"entrez_id"]})))
ENTREZID_case_earlyAD<-unique(unlist(lapply(Venn.case_earlyAD,function(x){ID_mapping[which(ID_mapping$GENE==x),"entrez_id"]})))

write.table(as.data.frame(ENTREZID_case_dementia),"./metadata/list_case_dementia_entrezid.txt",col.names = F,row.names = F,quote = F)
write.table(as.data.frame(ENTREZID_case_AD),"./metadata/list_case_AD_entrezid.txt",col.names = F,row.names = F,quote = F)
write.table(as.data.frame(ENTREZID_case_earlyAD),"./metadata/list_case_earlyAD_entrezid.txt",col.names = F,row.names = F,quote = F)

Venndata_dementia.caseonly_entrezid <- unlist(lapply(Venndata_dementia.caseonly$`3`, function(x){ID_mapping[which(ID_mapping$GENE %in% x),"entrez_id"]}))
Venndata_AD.caseonly_entrezid <- unique(unlist(lapply(Venndata_AD.caseonly$`3`, function(x){ID_mapping[which(ID_mapping$GENE %in% x),"entrez_id"]})))
Venndata_earlyAD.caseonly_entrezid <- unlist(lapply(Venndata_earlyAD.caseonly$`3`, function(x){ID_mapping[which(ID_mapping$GENE %in% x),"entrez_id"]}))


# Bar plot: mutation frequency -- cases in all groups

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


# ChowRuskey plot

## overlap of mutated genes and AD known genes

AD_knownsig_genes <- read.xlsx("./metadata/AD_knownsig_geneset_v2.xlsx")
ADknown_combined_allgroups <- list("AD known"=AD_knownsig_genes$Symbol,"dementia"=Venn.case_dementia,"AD"=Venn.case_AD,"earlyAD"=Venn.case_earlyAD)
V4a <- Venn(Sets = ADknown_combined_allgroups)
plot(V4a, type = "ChowRuskey", show = list(SetLabels = T))
partitions_ADknown_combined_allgroups <- get.venn.partitions(ADknown_combined_allgroups)
overlapped_gene_all <- partitions_ADknown_combined_allgroups[1,"..values.."]$`1`

## overlap of case only mutated genes and AD known genes

caseonly_combined_allgroups <- list("AD known"=AD_knownsig_genes$Entrez_ID,"dementia"=Venndata_dementia.caseonly_entrezid,"AD"=Venndata_AD.caseonly_entrezid,"earlyAD"=Venndata_earlyAD.caseonly_entrezid)
V4a_caseonly <- Venn(Sets = caseonly_combined_allgroups)
plot(V4a_caseonly, type = "ChowRuskey", show = list(SetLabels = T))
partitions_caseonly_combined_allgroups <- get.venn.partitions(caseonly_combined_allgroups)
overlapped_caseonly_all <- partitions_caseonly_combined_allgroups[1,"..values.."]$`1`

## overlap of burden genes and AD known genes

list_burden_dementia<-read.table("./metadata/burdengenes_dementia_list.txt")[,1]
list_burden_AD<-read.table("./metadata/burdengenes_AD_list.txt")[,1]
list_burden_earlyAD<-read.table("./metadata/burdengenes_ADearly_list.txt")[,1]
ADknown_burden_allgroups <- list("AD known"=AD_knownsig_genes$Entrez_ID,"dementia"=list_burden_dementia,"AD"=list_burden_AD,"earlyAD"=list_burden_earlyAD)
V4a_burden <- Venn(Sets = ADknown_burden_allgroups)
plot(V4a_burden, type = "ChowRuskey", show = list(SetLabels = T))
partitions_ADknown_burden_allgroups <- get.venn.partitions(ADknown_burden_allgroups)
overlapped_burden_all <- partitions_ADknown_burden_allgroups[1,"..values.."]$`1`

# save image

save.image("./VennPlot_TRAPDvariants_caseandcontrol.RData")

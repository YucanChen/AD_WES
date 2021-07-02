library(openxlsx)
GTEx_tissue_medianTPM <- read.table("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",sep = "\t",skip = 2,header = T)
burdengene_dementia <- read.csv("list_burdengenes_dementia.csv")
burdengene_AD <- read.csv("list_burdengenes_AD.csv")
burdengene_earlyAD <- read.csv("list_burdengenes_ADearly.csv")

# venn plot of the three burden genesets
library(Vennerable)
burdengene_allgroups <- list("dementia"=burdengene_dementia$GENE,"AD"=burdengene_AD$GENE,"earlyAD"=burdengene_earlyAD$GENE)
V3a <- Venn(Sets = burdengene_allgroups)
plot(V3a, type = "ChowRuskey", show = list(SetLabels = F))

# tissue-specific expression of genesets
expdf_ADknowngenes <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% AD_knownsig_genes$Symbol),];rownames(expdf_ADknowngenes) <- expdf_ADknowngenes$Description;expdf_ADknowngenes<-expdf_ADknowngenes[,3:ncol(expdf_ADknowngenes)]
expdf_burdengene_dementia <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% burdengene_dementia$GENE),];rownames(expdf_burdengene_dementia) <- expdf_burdengene_dementia$Description;expdf_burdengene_dementia<-expdf_burdengene_dementia[,3:ncol(expdf_burdengene_dementia)]
expdf_burdengene_AD <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% burdengene_AD$GENE),];rownames(expdf_burdengene_AD) <- expdf_burdengene_AD$Description;expdf_burdengene_AD<-expdf_burdengene_AD[,3:ncol(expdf_burdengene_AD)]
expdf_burdengene_earlyAD <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% burdengene_earlyAD$GENE),];rownames(expdf_burdengene_earlyAD) <- expdf_burdengene_earlyAD$Description;expdf_burdengene_earlyAD<-expdf_burdengene_earlyAD[,3:ncol(expdf_burdengene_earlyAD)]
expdf_caseonly_dementia <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% Venn.case_dementia),];rownames(expdf_caseonly_dementia) <- expdf_caseonly_dementia$Description;expdf_caseonly_dementia<-expdf_caseonly_dementia[,3:ncol(expdf_caseonly_dementia)]
expdf_caseonly_AD <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% Venn.case_AD),];rownames(expdf_caseonly_AD) <- expdf_caseonly_AD$Description;expdf_caseonly_AD<-expdf_caseonly_AD[,3:ncol(expdf_caseonly_AD)]
expdf_caseonly_earlyAD <- GTEx_tissue_medianTPM[which(GTEx_tissue_medianTPM$Description %in% Venn.case_earlyAD),];rownames(expdf_caseonly_earlyAD) <- expdf_caseonly_earlyAD$Description;expdf_caseonly_earlyAD<-expdf_caseonly_earlyAD[,3:ncol(expdf_caseonly_earlyAD)]

# heatmap plot
library(ComplexHeatmap)
ht_burdengene_dementia = Heatmap(as.matrix(expdf_burdengene_dementia), name = "burdengene dementia", cluster_columns = FALSE)
ht_burdengene_dementia
ht_burdengene_AD = Heatmap(as.matrix(expdf_burdengene_AD), name = "burdengene AD", cluster_columns = FALSE)
ht_burdengene_AD
ht_burdengene_earlyAD = Heatmap(as.matrix(expdf_burdengene_earlyAD), name = "burdengene earlyAD", cluster_columns = FALSE)
ht_burdengene_earlyAD
ht_ADknowngenes = Heatmap(as.matrix(expdf_ADknowngenes), name = "ADknowngenes", cluster_columns = FALSE)
ht_ADknowngenes

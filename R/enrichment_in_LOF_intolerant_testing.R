library(openxlsx)
# <Using variants from TRAPD outfiles>
# You need to load Venndata containing the case-only mutated genes before the analyses

# Part1: import pLI_data: pLI_B => all background; pLI_09 => pLI > 0.9
pLI_B <- read.csv("./gene_GeneticEvidence_MasterTable_foranalysis.csv")[,c(1:10,75,76)]
pLI_B <- pLI_B[complete.cases(pLI_B[,"ExAC.all.pLI...0.9"]), ]
pLI_09 <- pLI_B[which(pLI_B$ExAC.all.pLI...0.9==1),]

# Part2: <Enrichment analyses of case only mutated genes in intolerant gene list>
case_only_dementia_09 <- Venndata_dementia.caseonly$`3`[which(Venndata_dementia.caseonly$`3` %in% pLI_09$Symbol)]
case_only_AD_09 <- Venndata_AD.caseonly$`3`[which(Venndata_AD.caseonly$`3` %in% pLI_09$Symbol)]
case_only_earlyAD_09 <- Venndata_earlyAD.caseonly$`3`[which(Venndata_earlyAD.caseonly$`3` %in% pLI_09$Symbol)]
## Define the function
geneset_LOF_intolerant_testing <- function(geneset_character,pLI_09,pLI_B){
  case_09 <- geneset_character[(geneset_character %in% pLI_09$Symbol)]
  case_no09 <- geneset_character[(geneset_character %in% pLI_B$Symbol)]
  case_no09 <- case_no09[!(case_no09 %in% pLI_09$Symbol)]
  
  control_09_tmp <- pLI_B$Symbol[!(pLI_B$Symbol %in% geneset_character)]
  control_09 <- control_09_tmp[(control_09_tmp %in% pLI_09$Symbol)]
  control_no09 <- control_09_tmp[!(control_09_tmp %in% pLI_09$Symbol)]

  test.set <- c(length(case_09),length(control_09),length(case_no09),length(control_no09))

  test.dataframe <- matrix(test.set,nrow = 2,dimnames = list(c("case","control"),c("intolerant","tolerant")))
  test.fisher_result <- fisher.test(test.dataframe, alternative = "greater")
  test.fisher_result
  return(test.fisher_result)
}
## Caculate the enrichment of case only mutated genes in pLI_09 list
result_case_only_dementia<-geneset_LOF_intolerant_testing(case_only_dementia_09,pLI_09 = pLI_09,pLI_B = pLI_B)
result_case_only_AD<-geneset_LOF_intolerant_testing(case_only_AD_09,pLI_09 = pLI_09,pLI_B = pLI_B)
result_case_only_earlyAD<-geneset_LOF_intolerant_testing(case_only_earlyAD_09,pLI_09 = pLI_09,pLI_B = pLI_B)
## Check the significance
result_case_only_dementia
result_case_only_AD
result_case_only_earlyAD

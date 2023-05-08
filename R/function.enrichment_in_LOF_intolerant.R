geneset_LOF_intolerant_testing <- function(geneset_ENTREZID,pLI_09,pLI_B){
  case_09 <- geneset_ENTREZID[(geneset_ENTREZID %in% pLI_09$X.GeneID)]
  case_no09 <- geneset_ENTREZID[(geneset_ENTREZID %in% pLI_B$X.GeneID)]
  case_no09 <- case_no09[!(case_no09 %in% pLI_09$X.GeneID)]
  
  control_09_tmp <- pLI_B$X.GeneID[!(pLI_B$X.GeneID %in% geneset_ENTREZID)]
  control_09 <- control_09_tmp[(control_09_tmp %in% pLI_09$X.GeneID)]
  control_no09 <- control_09_tmp[!(control_09_tmp %in% pLI_09$X.GeneID)]

  test.set <- c(length(case_09),length(control_09),length(case_no09),length(control_no09))

  test.dataframe <- matrix(test.set,nrow = 2,dimnames = list(c("case","control"),c("intolerant","tolerant")))
  test.fisher_result <- fisher.test(test.dataframe, alternative = "greater")
  test.fisher_result
  return(test.fisher_result)
}

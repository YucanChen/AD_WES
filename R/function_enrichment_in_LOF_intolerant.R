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

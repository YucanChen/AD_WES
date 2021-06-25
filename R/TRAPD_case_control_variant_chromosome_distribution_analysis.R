library(openxlsx)

chr_case_dementia <- read.table("chr_case_dementia_outfile_snps.txt",sep = '\t',header = T)
chr_control_dementia <- read.table("chr_control_dementia_outfile_snps.txt",sep = '\t',header = T)
chr_case_AD <- read.table("chr_case_AD_outfile_snps.txt",sep = '\t',header = T)
chr_control_AD <- read.table("chr_control_AD_outfile_snps.txt",sep = '\t',header = T)
chr_case_earlyAD <- read.table("chr_case_earlyAD_outfile_snps.txt",sep = '\t',header = T)
chr_control_earlyAD <- read.table("chr_control_earlyAD_outfile_snps.txt",sep = '\t',header = T)

df_chr_dementia_test <- merge(chr_case_dementia,chr_control_dementia,by = "chr");colnames(df_chr_dementia_test) <- c("chr","case_snp","control_snp");rownames(df_chr_dementia_test) <- paste("chr",df_chr_dementia_test$chr,sep = "");df_chr_dementia_test <- df_chr_dementia_test[,-c(1)]
df_chr_AD_test <- merge(chr_case_AD,chr_control_AD,by = "chr");colnames(df_chr_AD_test) <- c("chr","case_snp","control_snp");rownames(df_chr_AD_test) <- paste("chr",df_chr_AD_test$chr,sep = "");df_chr_AD_test <- df_chr_AD_test[,-c(1)]
df_chr_earlyAD_test <- merge(chr_case_earlyAD,chr_control_earlyAD,by = "chr");colnames(df_chr_earlyAD_test) <- c("chr","case_snp","control_snp");rownames(df_chr_earlyAD_test) <- paste("chr",df_chr_earlyAD_test$chr,sep = "");df_chr_earlyAD_test <- df_chr_earlyAD_test[,-c(1)]

df_chr_dementia_test <- t(df_chr_dementia_test)
df_chr_AD_test <- t(df_chr_AD_test)
df_chr_earlyAD_test <- t(df_chr_earlyAD_test)

# chi-square test
chisq.test(df_chr_dementia_test)
chisq.test(df_chr_AD_test)
chisq.test(df_chr_earlyAD_test)

# individual chromosome level fisher test
chr_chisquare_test <- function(df){
  case_sum <- sum(df$case_snp)
  control_sum <- sum(df$control_snp)
  result_df <- data.frame()
  for(i in 1:length(df$chr)){
    test.set <- c(df[i,"case_snp"],df[i,"control_snp"],case_sum-df[i,"case_snp"],control_sum-df[i,"control_snp"])
    test.dataframe <- matrix(test.set,nrow = 2,dimnames = list(c("case","control"),c("chr_i","other_chrs")))
    result <- fisher.test(test.dataframe)
    chri_df <- data.frame(chr=df[i,"chr"],dds_ratio=as.numeric(result$estimate),p_value=result$p.value)
    result_df <- rbind(result_df,chri_df)
  }
  result_df$adj_pvalue <- p.adjust(result_df$p_value,method = "fdr")
  return(result_df)
}
result_dementia <- chr_chisquare_test(df_chr_dementia_test)
result_AD <- chr_chisquare_test(df_chr_AD_test)
result_earlyAD <- chr_chisquare_test(df_chr_earlyAD_test)
sheets = list("fisher_dementia"=result_dementia,"fisher_AD"=result_AD,"fisher_earlyAD"=result_earlyAD)
write.xlsx(sheets,"chromosome_distribution_fishertest_result_all.xlsx")

# output the dataframe for drawing the plots
df_chr_dementia_test <- as.data.frame(t(df_chr_dementia_test));df_chr_dementia_test$chr <- as.integer(unlist(lapply(rownames(df_chr_dementia_test), function(x){unlist(strsplit(x,"chr")[[1]][2])})))
df_chr_AD_test <- as.data.frame(t(df_chr_AD_test));df_chr_AD_test$chr <- as.integer(unlist(lapply(rownames(df_chr_AD_test), function(x){unlist(strsplit(x,"chr")[[1]][2])})))
df_chr_earlyAD_test <- as.data.frame(t(df_chr_earlyAD_test));df_chr_earlyAD_test$chr <- as.integer(unlist(lapply(rownames(df_chr_earlyAD_test), function(x){unlist(strsplit(x,"chr")[[1]][2])})))

write.csv(df_chr_dementia_test,"df_chr_dementia_test.csv")
write.csv(df_chr_AD_test,"df_chr_AD_test.csv")
write.csv(df_chr_earlyAD_test,"df_chr_earlyAD_test.csv")

# drawing forest plot
library("forestplot")
library(lattice)
## dementia
OR_interval_df_dementia <- 
  structure(list(
    mean  = c(NA,result_dementia$odds_ratio), 
    lower = c(NA,result_dementia$conf.int_lower),
    upper = c(NA,result_dementia$conf.int_upper)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -23L),
    class = "data.frame")
tabletext_dementia<-cbind(
  c("Chromosome",paste("chr",result_dementia$chr,sep = "")),
  c("OR",sprintf("%0.2f",result_dementia$odds_ratio)),
  c("CI",paste("(",paste(sprintf("%0.2f",result_dementia$conf.int_lower),sprintf("%0.2f",result_dementia$conf.int_upper),sep = "-"),")",sep = "")),
  c("p.value",sprintf("%0.3f",result_dementia$p_value)),
  c("p.adjust",sprintf("%0.3f",result_dementia$adj_pvalue)))
trellis.device(device = "windows", height = 25, width = 40, color = TRUE)
forestplot(tabletext_dementia,
           graph.pos = 2,
           hrzl_lines = list("2" = gpar(lty=2),
                            "24" = gpar(lwd=1, lty=2, col = "#000044")),
           OR_interval_df_dementia,new_page = TRUE,
           clip=c(0.1,2.5), 
           xlog=TRUE,
           col=fpColors(box="royalblue",line="darkblue", hrz_lines = "#444444")
          )
## AD
OR_interval_df_AD <- 
  structure(list(
    mean  = c(NA,result_AD$odds_ratio), 
    lower = c(NA,result_AD$conf.int_lower),
    upper = c(NA,result_AD$conf.int_upper)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -23L),
    class = "data.frame")
tabletext_AD<-cbind(
  c("Chromosome",paste("chr",result_AD$chr,sep = "")),
  c("OR",sprintf("%0.2f",result_AD$odds_ratio)),
  c("CI",paste("(",paste(sprintf("%0.2f",result_AD$conf.int_lower),sprintf("%0.2f",result_AD$conf.int_upper),sep = "-"),")",sep = "")),
  c("p.value",sprintf("%0.3f",result_AD$p_value)),
  c("p.adjust",sprintf("%0.3f",result_AD$adj_pvalue)))
trellis.device(device = "windows", height = 25, width = 40, color = TRUE)
forestplot(tabletext_AD,
           graph.pos = 2,
           hrzl_lines = list("2" = gpar(lty=2),
                             "24" = gpar(lwd=1, lty=2, col = "#000044")),
           OR_interval_df_AD,new_page = TRUE,
           clip=c(0.1,2.5), 
           xlog=TRUE,
           col=fpColors(box="royalblue",line="darkblue", hrz_lines = "#444444")
)
## AD early
OR_interval_df_earlyAD <- 
  structure(list(
    mean  = c(NA,result_earlyAD$odds_ratio), 
    lower = c(NA,result_earlyAD$conf.int_lower),
    upper = c(NA,result_earlyAD$conf.int_upper)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -23L),
    class = "data.frame")
tabletext_earlyAD<-cbind(
  c("Chromosome",paste("chr",result_earlyAD$chr,sep = "")),
  c("OR",sprintf("%0.2f",result_earlyAD$odds_ratio)),
  c("CI",paste("(",paste(sprintf("%0.2f",result_earlyAD$conf.int_lower),sprintf("%0.2f",result_earlyAD$conf.int_upper),sep = "-"),")",sep = "")),
  c("p.value",sprintf("%0.3f",result_earlyAD$p_value)),
  c("p.adjust",sprintf("%0.3f",result_earlyAD$adj_pvalue)))
trellis.device(device = "windows", height = 25, width = 40, color = TRUE)
forestplot(tabletext_earlyAD,
           graph.pos = 2,
           hrzl_lines = list("2" = gpar(lty=2),
                             #"11" = gpar(lwd=1, columns=c(1:3,5), col = "#000044"),
                             "24" = gpar(lwd=1, lty=2, col = "#000044")),
           OR_interval_df_earlyAD,new_page = TRUE,
           clip=c(0.1,2.5), 
           xlog=TRUE,
           col=fpColors(box="royalblue",line="darkblue", hrz_lines = "#444444")
)

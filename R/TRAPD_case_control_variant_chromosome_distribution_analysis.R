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

chisq.test(df_chr_dementia_test)
chisq.test(df_chr_AD_test)
chisq.test(df_chr_earlyAD_test)

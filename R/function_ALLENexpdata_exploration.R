#' @import dplyr,magrittr,tidyverse,limma,sva

# DE analysis-expression matrix preparation

get_exprmat_ALLEN<-function(exp_matrix,exp_colsamples,brain_region,disease_type,geneset_entrezid="NA"){
  # prepare expression matrix
  
  ## filter brain region
  if(disease_type=="AD"){
    exp_colsample_brainregion<-exp_colsamples %>% filter(structure_acronym==brain_region) %>% filter(dsm_iv_clinical_diagnosis=="Alzheimers Disease Type"|dsm_iv_clinical_diagnosis=="No Dementia")
  }else if(disease_type=="dementia"){
    exp_colsample_brainregion<-exp_colsamples %>% filter(structure_acronym==brain_region)
  }
  exp_colsample_brainregion$rnaseq_profile_id<-paste("X",exp_colsample_brainregion$rnaseq_profile_id,sep = '')
  exp_matrix %<>% select(all_of(exp_colsample_brainregion$rnaseq_profile_id))
  ## filter low expressed genes
  exp_matrix %<>% cbind(.,sum=rowSums(.))
  exp_matrix %<>% rownames_to_column() %>% filter(sum > 1) %>% select(-c("sum")) %<>% column_to_rownames()
  
  # perform DE analysis
  if(disease_type=="AD"){
    exp_condition<-unlist(lapply(exp_colsample_brainregion$nincds_arda_diagnosis, function(x){if(x=="No Dementia"){return("NoDementia")}else{return("AlzheimersDisease")}}))
  }else if(disease_type=="dementia"){
    exp_condition=unlist(lapply(exp_colsample_brainregion$nincds_arda_diagnosis, function(x){if(x=="No Dementia"){return("Normal")}else{return("Dementia")}}))
  }
  condition <- factor(exp_condition)
  pheno=data.frame("sample"=colnames(exp_matrix),"condition"=condition,"batch"=exp_colsample_brainregion$structure_hemisphere)
  mod = model.matrix(~as.factor(condition), data=pheno)
  mod0 = model.matrix(~as.factor(batch),data=pheno)
  exp_matrix %<>% as.matrix
  n.sv = num.sv(exp_matrix,mod,method="leek")
  print(n.sv)
  design <- model.matrix(~0+condition, data = pheno)
  if(n.sv>0){
    svseq <- sva(exp_matrix, mod, n.sv = n.sv)
    combat_edata <- removeBatchEffect(exp_matrix, covariates=svseq$sv[,1],design = design)
    dist_mat_combat <- dist(t(combat_edata))
    clustering_combat <- hclust(dist_mat_combat, method = "ward.D")
    plot(clustering_combat, labels = pheno$condition)
  }else{combat_edata <- exp_matrix}
  ## Note: the matrix has been log2 preprocessed, so svaseq is not needed. Moreover, mod0 is conflict with sv in batch, so mod0 is not used in this situation
  
  # detect batch effect (visualization)
  
  ## check the original data
  dist_mat <- dist(t(exp_matrix))
  clustering <- hclust(dist_mat, method = "ward.D")
  plot(clustering, labels = pheno$condition)
  colramp = colorRampPalette(c(3,"white",2))(length(combat_edata[1,]))
  plot(density(exp_matrix[,1]),col=colramp[1],lwd=3,ylim=c(0,.80))
  for(i in 2:length(exp_matrix[1,])){lines(density(exp_matrix[,i]),lwd=3,col=colramp[i])}
  
  ## delete batch effect: Hemisphere
  batch = pheno$batch
  combat_edata <- removeBatchEffect(combat_edata,batch = pheno$batch,design = design)
  dist_mat_combat <- dist(t(combat_edata))
  clustering_combat <- hclust(dist_mat_combat, method = "ward.D")
  plot(clustering_combat, labels = pheno$condition)
  
  ## delete batch effect: age
  pheno=data.frame("sample"=colnames(combat_edata),"condition"=condition,"batch"=as.factor(exp_colsample_brainregion$age))
  batch = pheno$batch
  combat_edata <- removeBatchEffect(combat_edata,batch = pheno$batch,design = design)
  dist_mat_combat <- dist(t(combat_edata))
  clustering_combat <- hclust(dist_mat_combat, method = "ward.D")
  plot(clustering_combat, labels = pheno$condition)
  
  ## visualize for rechecking
  colramp = colorRampPalette(c(3,"white",2))(length(combat_edata[1,]))
  plot(density(combat_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.80))
  for(i in 2:length(combat_edata[1,])){lines(density(combat_edata[,i]),lwd=3,col=colramp[i])}
  
  print("Finished")
  
  combat_edata %<>% as.data.frame
  if(geneset_entrezid=="NA"){
    return(combat_edata %<>% rbind(.,nincds_arda_diagnosis=exp_condition))
  }else{
    combat_edata.filtered<-combat_edata[which(rownames(combat_edata) %in% geneset_entrezid),]
    return(combat_edata.filtered %<>% rbind(.,nincds_arda_diagnosis=exp_condition))
  }
}


# obtain significantly DE genes

get_signifgenes<-function(brain_region_combatedata,brain_region,disease_type,fdr_threshold){
  # prepare expression matrix
  exp_condition<-brain_region_combatedata[which(rownames(brain_region_combatedata)=="nincds_arda_diagnosis"),] %>% as.character
  zcore_combat_exp<-brain_region_combatedata[which(rownames(brain_region_combatedata)!="nincds_arda_diagnosis"),]
  rownames(zcore_combat_exp)=unlist(lapply(rownames(zcore_combat_exp), function(x){list_caseonly_and_burden_all[which(list_caseonly_and_burden_all$ENTREZID==x),"Symbol"]}))
  ## Note: list_caseonly_and_burden_all: here is the table containing all testing genes for ID conversion
  zcore_combat_exp<-zcore_combat_exp[apply(zcore_combat_exp,1,var)!=0,]
  
  # scaling for Z-score
  zcore_combat_exp.row<-rownames(zcore_combat_exp)
  zcore_combat_exp %<>% lapply(.,as.numeric) %>% as.data.frame %>% apply(.,1,scale) %>% Matrix::t() %>% as.data.frame
  rownames(zcore_combat_exp)<-zcore_combat_exp.row
  colnames(zcore_combat_exp)<-exp_condition
  
  # DE analysis
  zcore_combat_exp.pvalue=c()
  if(disease_type=="AD"){
    for(i in 1:length(zcore_combat_exp[,1])){
      A=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="AlzheimersDisease")])
      B=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="NoDementia")])
      result=t.test(A,B,alternative = c("two.sided"))
      zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
    }
  }else if(disease_type=="dementia"){
    for(i in 1:length(zcore_combat_exp[,1])){
      A=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="Dementia")])
      B=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="Normal")])
      result=t.test(A,B,alternative = c("two.sided"))
      zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
    }
  }
  signif_genes=cbind(zcore_combat_exp,as.data.frame(zcore_combat_exp.pvalue) %>% dplyr::rename(c(p.value=zcore_combat_exp.pvalue)))
  signif_genes$padjust_fdr=p.adjust(signif_genes$p.value,method = "fdr")
  signif_genes<-signif_genes[order(signif_genes$padjust_fdr),]
  signif_genes$sig_DE<-lapply(signif_genes$padjust_fdr,function(x){if(x < fdr_threshold){"TRUE"}else{"FALSE"}}) %>% unlist
  
  if("TRUE" %in% signif_genes$sig_DE){
    print("exist!")
    return(signif_genes)
  }else{
    print("no significantly DE genes")
  }
}


# plot

#' @import circlize,ComplexHeatmap

## heatmap plot

draw_heatmap_zscore<-function(zcore_signifgenes,disease_type,brain_region,fdr_threshold="NA"){
  col_fun = colorRamp2(c(-4,-2,0,2,4),c("#330033","#660066", "#CC0000", "yellow","#FFFF99"))
  
  if(fdr_threshold=="NA"){
    matrix_heatplot<-as.matrix(zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)])
  }else{
    matrix_heatplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)] %>% cbind(.,padjust_fdr=zcore_signifgenes$padjust_fdr) %>% rownames_to_column() %>% filter(padjust_fdr < fdr_threshold) %>% column_to_rownames() %>% select(-c("padjust_fdr")) %>% as.matrix
  }
  if(disease_type=="AD"){
    ha = HeatmapAnnotation(type = unlist(lapply(colnames(matrix_heatplot), function(x){strsplit(x,split = '\\.')[[1]][1]})),
                           col = list(type = c("AlzheimersDisease" = "#8B3A62", "NoDementia" = "#B0C4DE")),
                           annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(matrix_heatplot)) + unit(2, "mm")))
  }else if(disease_type=="dementia"){
    ha = HeatmapAnnotation(type = unlist(lapply(colnames(matrix_heatplot), function(x){strsplit(x,split = '\\.')[[1]][1]})),
                           col = list(type = c("Dementia" = "#8B3A62", "Normal" = "#B0C4DE")),
                           annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(matrix_heatplot)) + unit(2, "mm")))
  }
  
  Heatmap(matrix_heatplot, show_column_names = TRUE, bottom_annotation = ha,col=col_fun,name = "expression",column_title = paste(disease_type,brain_region,sep = " vs Norm: "))
}

## violin plot

draw_violin_plots<-function(zcore_signifgenes,disease_type,brain_region,fdr_threshold="NA"){
  if(fdr_threshold=="NA"){
    matrix_violinplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)]
  }else{
    matrix_violinplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)] %>% cbind(.,padjust_fdr=zcore_signifgenes$padjust_fdr) %>% rownames_to_column() %>% filter(padjust_fdr < fdr_threshold) %>% column_to_rownames() %>% select(-c("padjust_fdr"))
  }
  matrix_violinplot$gene<-rownames(matrix_violinplot)
  matrix_violinplot %<>% melt
  colnames(matrix_violinplot)<-c("gene","group","zscore_expression")
  matrix_violinplot$group %<>% as.character %>% lapply(.,function(x){strsplit(x,split = '\\.')[[1]][1]}) %>% unlist
  
  ggviolin(matrix_violinplot,x="gene",y="zscore_expression",
           color = "group",
           palette = c("#8B3A62","#4F94CD"),
           alpha=0.7,
           linetype=1,
           trim=F,
           size=1,
           width=1, #violin width
           draw_quantiles = 0.5, #median value
           add = "boxplot",
           add.params = list(shape="group",color="group",alpha=0.6),
           error.plot = "errorbar",
           position=position_dodge(1),
           xlab="burden genes",
           ylab="expression z-scores",
           title=paste(disease_type,brain_region,sep = " vs Norm: "),
           ggtheme=theme_pubr(),
  )
}

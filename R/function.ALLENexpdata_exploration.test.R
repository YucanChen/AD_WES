#' @import dplyr,magrittr,tidyverse,limma,sva

# DE analysis-expression matrix preparation

get_exprmat_ALLEN<-function(exp_matrix,exp_colsamples,brain_region,disease_type,geneset_entrezid="NA"){
  # prepare expression matrix
  
  ## filter brain region
  if(disease_type=="AD"){
    exp_colsample_brainregion<-exp_colsamples %>% dplyr::filter(structure_acronym==brain_region) %>% dplyr::filter(dsm_iv_clinical_diagnosis=="Alzheimers Disease Type"|dsm_iv_clinical_diagnosis=="No Dementia")
  }else if(disease_type=="dementia"){
    exp_colsample_brainregion<-exp_colsamples %>% dplyr::filter(structure_acronym==brain_region)
  }
  exp_colsample_brainregion$rnaseq_profile_id<-paste("X",exp_colsample_brainregion$rnaseq_profile_id,sep = '')
  exp_matrix %<>% dplyr::select(all_of(exp_colsample_brainregion$rnaseq_profile_id))
  ## filter low expressed genes
  exp_matrix %<>% cbind(.,sum=rowSums(.))
  exp_matrix %<>% rownames_to_column()
  exp_matrix %<>% dplyr::filter(sum > 1) %>% dplyr::select(-c("sum"))
  exp_matrix %<>% column_to_rownames()
  
  # perform batch effect correction
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
    for(i in 1:n.sv){
      combat_edata <- removeBatchEffect(exp_matrix, covariates=svseq$sv[,i],design = design)
    }
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

get_signifgenes<-function(brain_region_combatedata,brain_region,disease_type,DE_method="t.test",fdr_threshold,list_mut_all.Sym_ENTREZID){
  # prepare expression matrix
  exp_condition<-brain_region_combatedata[which(rownames(brain_region_combatedata)=="nincds_arda_diagnosis"),] %>% as.character
  zcore_combat_exp<-brain_region_combatedata[which(rownames(brain_region_combatedata)!="nincds_arda_diagnosis"),]
  rownames(zcore_combat_exp)=unlist(lapply(rownames(zcore_combat_exp), function(x){list_mut_all.Sym_ENTREZID[which(list_mut_all.Sym_ENTREZID$ENTREZID==x),"GENE"]}))
  
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
    for(i in 1:nrow(zcore_combat_exp)){
      A=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="AlzheimersDisease")])
      B=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="NoDementia")])
      if(DE_method == 't.test'){
        result=t.test(A,B,alternative = c("two.sided"))
        zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
      }
      if(DE_method == 'wilcox.test'){
        result=wilcox.test(A,B,alternative = c("two.sided"))
        zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
      }
    }
  }else if(disease_type=="dementia"){
    for(i in 1:nrow(zcore_combat_exp)){
      A=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="Dementia")])
      B=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="Normal")])
      if(DE_method == 't.test'){
        result=t.test(A,B,alternative = c("two.sided"))
        zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
      }
      if(DE_method == 'wilcox.test'){
        result=wilcox.test(A,B,alternative = c("two.sided"))
        zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
      }
    }
  }
  signif_genes=cbind(zcore_combat_exp,as.data.frame(zcore_combat_exp.pvalue) %>% dplyr::rename(c(p.value=zcore_combat_exp.pvalue)))
  signif_genes$padjust_fdr=p.adjust(signif_genes$p.value,method = "fdr")
  signif_genes<-signif_genes[order(signif_genes$padjust_fdr),]
  signif_genes$sig_DE<-lapply(signif_genes$padjust_fdr,function(x){if(x < fdr_threshold){"TRUE"}else{"FALSE"}}) %>% unlist
  rownames(signif_genes)<-zcore_combat_exp.row
  
  if("TRUE" %in% signif_genes$sig_DE){
    print("exist!")
    return(signif_genes)
  }else{
    print("no significantly DE genes")
  }
}


get_signifgenes.test<-function(brain_region_combatedata,brain_region,disease_type,DE_method="t.test",fdr_threshold){
  # prepare expression matrix
  exp_condition<-brain_region_combatedata[which(rownames(brain_region_combatedata)=="nincds_arda_diagnosis"),] %>% as.character
  zcore_combat_exp<-brain_region_combatedata[which(rownames(brain_region_combatedata)!="nincds_arda_diagnosis"),]
  ## Note: list_caseonly_and_burden_all: here is the table containing all testing genes for ID conversion
  zcore_combat_exp<-zcore_combat_exp[apply(zcore_combat_exp,1,var)!=0,]
  zcore_combat_exp.row<-rownames(zcore_combat_exp)
  
  # DE analysis
  zcore_combat_exp.pvalue=c()
  if(disease_type=="AD"){
    if((DE_method == 't.test')|(DE_method == 'wilcox.test')){
      # scaling for Z-score
      zcore_combat_exp %<>% lapply(.,as.numeric) %>% as.data.frame %>% apply(.,1,scale) %>% Matrix::t() %>% as.data.frame
      rownames(zcore_combat_exp)<-zcore_combat_exp.row
      colnames(zcore_combat_exp)<-exp_condition
      
      for(i in 1:nrow(zcore_combat_exp)){
        A=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="AlzheimersDisease")])
        B=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="NoDementia")])
        if(DE_method == 't.test'){
          result=t.test(A,B,alternative = c("two.sided"))
          zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
        }else{
          result=wilcox.test(A,B,alternative = c("two.sided"))
          zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
        }
      }
      signif_genes=cbind(zcore_combat_exp,as.data.frame(zcore_combat_exp.pvalue) %>% dplyr::rename(c(p.value=zcore_combat_exp.pvalue)))
      signif_genes$padjust_fdr=p.adjust(signif_genes$p.value,method = "fdr")
      rownames(signif_genes)<-zcore_combat_exp.row
      signif_genes<-signif_genes[order(signif_genes$padjust_fdr),]
      signif_genes$sig_DE<-lapply(signif_genes$padjust_fdr,function(x){if(x < fdr_threshold){"TRUE"}else{"FALSE"}}) %>% unlist
    }
    
    if(DE_method == 'limma'){
      u<-unique(exp_condition)
      f<-factor(exp_condition,levels = u)
      design<-model.matrix(~0+f)
      colnames(design)<-u
      cont.matrix<-makeContrasts(AlzheimersDiseasevsNoDementia=AlzheimersDisease-NoDementia, levels=design)
      
      colnames(zcore_combat_exp)<-paste(exp_condition,seq(1:length(exp_condition)),sep = '.')
      eset<-ExpressionSet(assayData=(zcore_combat_exp %>% apply(.,2,as.numeric)),phenoData = annotatedDataFrameFrom(zcore_combat_exp %>% as.matrix(), byrow=FALSE),featureData = annotatedDataFrameFrom(zcore_combat_exp %>% as.matrix(), byrow=TRUE))
      fit<-lmFit(eset,design)
      fit2<-contrasts.fit(fit, cont.matrix)
      fit2<-eBayes(fit2)
      top<-topTable(fit2, adjust="fdr",number=Inf,sort.by="P")
      print(paste('limma DE gene number:',sum(top$adj.P.Val<0.05),sep = ' '))
      
      signif_genes<-top %>% dplyr::rename(.,padjust_fdr=adj.P.Val,p.value=P.Value)
      signif_genes<-signif_genes[order(signif_genes$padjust_fdr),]
      signif_genes$sig_DE<-lapply(signif_genes$padjust_fdr,function(x){if(x < fdr_threshold){"TRUE"}else{"FALSE"}}) %>% unlist
    }
    
    #if(DE_method == 'lefse'){
    #  assaysData<-zcore_combat_exp %>% apply(.,2,as.numeric)
    #  colnames(assaysData)<-paste(exp_condition,seq(1:length(exp_condition)),sep = '.')
    #  colData<-DataFrame(disease_status=exp_condition,row.names=paste(exp_condition,seq(1:length(exp_condition)),sep = '.'))
    #  se0<-SummarizedExperiment(assays=SimpleList(expr=assaysData),
    #                            colData=colData)
    #  res_group.lefse<-lefser(se0,groupCol = 'disease_status',kruskal.threshold = 0.1,wilcox.threshold = 0.1,lda.threshold = 1.5)
      
    #  plt.lefse_1<-(se0@colData$disease_status %>% unique() %>% sort())[2]
    #  plt.lefse_0<-(se0@colData$disease_status %>% unique() %>% sort())[1]
    #  group<-ifelse(res_group.lefse$scores > 0, plt.lefse_1, plt.lefse_0)
    #  res_group.lefse$group<-factor(group,levels = c(plt.lefse_0,plt.lefse_1))
    #}
    
  }else if(disease_type=="dementia"){
    
    if((DE_method == 't.test')|(DE_method == 'wilcox.test')){
      # scaling for Z-score
      zcore_combat_exp %<>% lapply(.,as.numeric) %>% as.data.frame %>% apply(.,1,scale) %>% Matrix::t() %>% as.data.frame
      rownames(zcore_combat_exp)<-zcore_combat_exp.row
      colnames(zcore_combat_exp)<-exp_condition
      
      for(i in 1:nrow(zcore_combat_exp)){
        A=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="Dementia")])
        B=as.numeric(zcore_combat_exp[i,which(colnames(zcore_combat_exp)=="Normal")])
        if(DE_method == 't.test'){
          result=t.test(A,B,alternative = c("two.sided"))
          zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
        }else{
          result=wilcox.test(A,B,alternative = c("two.sided"))
          zcore_combat_exp.pvalue=c(zcore_combat_exp.pvalue,result$p.value)
        }
      }
      signif_genes=cbind(zcore_combat_exp,as.data.frame(zcore_combat_exp.pvalue) %>% dplyr::rename(c(p.value=zcore_combat_exp.pvalue)))
      signif_genes$padjust_fdr=p.adjust(signif_genes$p.value,method = "fdr")
      rownames(signif_genes)<-zcore_combat_exp.row
      signif_genes<-signif_genes[order(signif_genes$padjust_fdr),]
      signif_genes$sig_DE<-lapply(signif_genes$padjust_fdr,function(x){if(x < fdr_threshold){"TRUE"}else{"FALSE"}}) %>% unlist
    }
      
    if(DE_method == 'limma'){
      u<-unique(exp_condition)
      f<-factor(exp_condition,levels = u)
      design<-model.matrix(~0+f)
      colnames(design)<-u
      cont.matrix<-makeContrasts(DementiavsNormal=Dementia-Normal, levels=design)
      
      colnames(zcore_combat_exp)<-paste(exp_condition,seq(1:length(exp_condition)),sep = '.')
      eset<-ExpressionSet(assayData=(zcore_combat_exp %>% apply(.,2,as.numeric)),phenoData = annotatedDataFrameFrom(zcore_combat_exp %>% as.matrix(), byrow=FALSE),featureData = annotatedDataFrameFrom(zcore_combat_exp %>% as.matrix(), byrow=TRUE))
      fit<-lmFit(eset,design)
      fit2<-contrasts.fit(fit, cont.matrix)
      fit2<-eBayes(fit2)
      top<-topTable(fit2, adjust="fdr",number=Inf,sort.by="P")
      print(paste('limma DE gene number:',sum(top$adj.P.Val<0.05),sep = ' '))
      
      signif_genes<-top %>% dplyr::rename(.,padjust_fdr=adj.P.Val,p.value=P.Value)
      signif_genes<-signif_genes[order(signif_genes$padjust_fdr),]
      signif_genes$sig_DE<-lapply(signif_genes$padjust_fdr,function(x){if(x < fdr_threshold){"TRUE"}else{"FALSE"}}) %>% unlist
    }
  }
  
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
    matrix_heatplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)] %>% cbind(.,padjust_fdr=zcore_signifgenes$padjust_fdr) %>% rownames_to_column() %>% dplyr::filter(padjust_fdr < fdr_threshold) %>% column_to_rownames() %>% dplyr::select(-c("padjust_fdr")) %>% as.matrix
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
    matrix_violinplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)] %>% cbind(.,padjust_fdr=zcore_signifgenes$padjust_fdr) %>% rownames_to_column() %>% dplyr::filter(padjust_fdr < fdr_threshold) %>% column_to_rownames() %>% dplyr::select(-c("padjust_fdr"))
  }
  matrix_violinplot$gene<-rownames(matrix_violinplot)
  matrix_violinplot %<>% melt
  colnames(matrix_violinplot)<-c("gene","group","zscore_expression")
  matrix_violinplot$group %<>% as.character %>% lapply(.,function(x){strsplit(x,split = '\\.')[[1]][1]}) %>% unlist
  #注意这里x,y要用 "" string来表示
  
  ggviolin(matrix_violinplot,x="gene",y="zscore_expression",
           color = "group",
           palette = c("#8B3A62","#4F94CD"),
           alpha=0.7, #透明度
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
  )#+stat_compare_means(aes(group=group), label = "p.signif",method = "t.test")
}




# plot (test)

#' @import circlize,ComplexHeatmap

## heatmap plot

draw_heatmap_zscore.test<-function(zcore_signifgenes,disease_type,brain_region,fdr_threshold="NA",gene_list.entriz="NA"){
  col_fun = colorRamp2(c(-4,-2,0,2,4),c("#330033","#660066", "#CC0000", "yellow","#FFFF99"))
  zcore_signifgenes.row<-zcore_signifgenes %>% rownames()
  
  if(gene_list.entriz!="NA"){
    zcore_signifgenes %<>% tibble::add_column(ENTREZID=zcore_signifgenes.row) %>% dplyr::filter(ENTREZID %in% gene_list.entriz$ENTREZID)
    zcore_signifgenes.entriz<-zcore_signifgenes$ENTREZID
    rownames(zcore_signifgenes)<-lapply(zcore_signifgenes.entriz,function(x){gene_list.entriz[which(gene_list.entriz$ENTREZID==x),'GENE']}) %>% unlist
    zcore_signifgenes %<>% dplyr::select(!("ENTREZID"))
  }
  
  if(fdr_threshold=="NA"){
    matrix_heatplot<-as.matrix(zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)])
  }else{
    matrix_heatplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)] %>% cbind(.,padjust_fdr=zcore_signifgenes$padjust_fdr) %>% rownames_to_column() %>% dplyr::filter(padjust_fdr < fdr_threshold) %>% column_to_rownames() %>% dplyr::select(-c("padjust_fdr")) %>% as.matrix
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

draw_violin_plots.test<-function(zcore_signifgenes,disease_type,brain_region,fdr_threshold="NA",gene_list.entriz="NA"){
  zcore_signifgenes.row<-zcore_signifgenes %>% rownames()
  
  if(gene_list.entriz!="NA"){
    zcore_signifgenes %<>% tibble::add_column(ENTREZID=zcore_signifgenes.row) %>% dplyr::filter(ENTREZID %in% gene_list.entriz$ENTREZID)
    zcore_signifgenes.entriz<-zcore_signifgenes$ENTREZID
    rownames(zcore_signifgenes)<-lapply(zcore_signifgenes.entriz,function(x){gene_list.entriz[which(gene_list.entriz$ENTREZID==x),'GENE']}) %>% unlist
    zcore_signifgenes %<>% dplyr::select(!("ENTREZID"))
  }
  
  if(fdr_threshold=="NA"){
    matrix_violinplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)]
  }else{
    matrix_violinplot<-zcore_signifgenes[,1:(ncol(zcore_signifgenes)-3)] %>% cbind(.,padjust_fdr=zcore_signifgenes$padjust_fdr) %>% rownames_to_column() %>% dplyr::filter(padjust_fdr < fdr_threshold) %>% column_to_rownames() %>% dplyr::select(-c("padjust_fdr"))
  }
  matrix_violinplot$gene<-rownames(matrix_violinplot)
  matrix_violinplot %<>% melt
  colnames(matrix_violinplot)<-c("gene","group","zscore_expression")
  matrix_violinplot$group %<>% as.character %>% lapply(.,function(x){strsplit(x,split = '\\.')[[1]][1]}) %>% unlist
  matrix_violinplot$padjust_fdr<-lapply(matrix_violinplot$gene,function(x){zcore_signifgenes[which(rownames(zcore_signifgenes)==x),'padjust_fdr']}) %>% unlist
  #注意这里x,y要用 "" string来表示
  
  p.violin<- 
  ggviolin(matrix_violinplot,x="gene",y="zscore_expression",
           color = "group",
           palette = c("#8B3A62","#4F94CD"),
           alpha=0.7, #透明度
           linetype=1,
           trim=F,
           size=1,
           width=1, #violin width
           draw_quantiles = 0.5, #median value
           add = "boxplot",
           add.params = list(shape="group",color="group",alpha=0.6),
           error.plot = "errorbar",
           position=position_dodge(1),
           xlab="",
           ylab="expression z-scores",
           title=paste(disease_type,brain_region,sep = " vs Norm: "),
           ggtheme=theme_pubr(),
  ) + theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0.05))
  
  if(fdr_threshold=="NA"){
    print(p.violin)
  }else{
    label.padjust_fdr<-data.frame(gene = matrix_violinplot$gene %>% unique(), zscore_expression = rep(matrix_violinplot$zscore_expression %>% max() %>% add(.,0.6),matrix_violinplot$gene %>% unique() %>% length()), padjust_fdr = matrix_violinplot$padjust_fdr %>% unique())
    label.asterisk<-round(label.padjust_fdr$padjust_fdr,digits = 3)
    #if(fdr_threshold<=0.001){
    #  label.asterisk<-rep('***',nrow(label.padjust_fdr))
    #}
    #if((fdr_threshold<=0.05)&(fdr_threshold>0.001)){
    #  label.asterisk<-lapply(label.padjust_fdr$padjust_fdr, function(x){if(x<0.001){'***'}else if((x>=0.001)&(x<0.05)){'**'}}) %>% unlist
    #}
    #if((fdr_threshold<=0.1)&(fdr_threshold>0.05)){
    #  label.asterisk<-lapply(label.padjust_fdr$padjust_fdr, function(x){if(x<0.001){'***'}else if((x>=0.001)&(x<0.05)){'**'}else if((x>=0.05)&(x<0.1)){'*'}}) %>% unlist
    #}
    #if(fdr_threshold>0.1){
    #  label.asterisk<-lapply(label.padjust_fdr$padjust_fdr, function(x){if(x<0.001){'***'}else if((x>=0.001)&(x<0.05)){'**'}else if((x>=0.05)&(x<0.1)){'*'}else if(x>=0.1){''}}) %>% unlist
    #}
    p.violin<-p.violin + geom_text(data = label.padjust_fdr, label = label.asterisk,colour='dodgerblue4',size=3.5)
    print(p.violin)
  }
}

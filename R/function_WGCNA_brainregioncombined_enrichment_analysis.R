get_WGCNAconsensusmodule<-function(burdenresult,geneinfo){
  #adding the title using the celltype marker list
  #geneinfo must have been deleted with grey module
  module_names<-unique(geneinfo$ModuleColor)
  module_df<-burdenresult$GENE
  #build 0 matrix using the matrix
  for(i in 1:length(module_names)){
    assign(module_names[i],rep(0,length(burdenresult[,1])))
    module_df<-cbind(module_df,get(module_names[i]))
  }
  colnames(module_df)<-c("GENE",module_names)
  burdenresult<-merge(burdenresult, module_df, by = "GENE")
  return(burdenresult)
}

get_module_annotation<-function(burdenresult,geneinfo){
  module_list=unique(geneinfo$ModuleColor)
  for(i in 1:length(module_list)){
    module_df<-geneinfo[which(geneinfo$ModuleColor==module_list[i]),]
    burdenresult[,19+i]<-unlist(lapply(burdenresult$entrezgene_id,function(x){if(x %in% module_df$EntrezID){1}else{0}}))
  }
  return(burdenresult)
}

get_moduleEnrichment<-function(burdenresult){
  all_fit_disease<-data.frame()
  burdenresult_colnames<-colnames(burdenresult)
  for(i in 20:length(burdenresult[1,])){
    burdenresult.fit<-glm(disease~burdenresult[,i]+percentage_gene_gc_content+transcript_length,data = burdenresult)
    tmp<-summary(burdenresult.fit)
    tmp<-as.data.frame(tmp$coefficients)
    tmp<-tmp["burdenresult[, i]",]
    rownames(tmp)<-burdenresult_colnames[i]
    all_fit_disease<-rbind(all_fit_disease,tmp)
  }
  all_fit_disease$modulecolor<-rownames(all_fit_disease)
  all_fit_disease$p_adjust<-p.adjust(all_fit_disease$`Pr(>|t|)`,method = "fdr")
  return(all_fit_disease)
}

get_diseasestate<-function(burdenresult,P_DOM_adjustFDR,CASE_TOTAL_AC){
  disease<-c()
  for(i in 1:length(burdenresult$entrezgene_id)){
    if((burdenresult[i,"P_DOM_adjustFDR"]<P_DOM_adjustFDR)|(burdenresult[i,"CASE_TOTAL_AC"]>CASE_TOTAL_AC)){
      disease<-c(disease,1)
    }else{disease<-c(disease,0)}
  }
  return(disease)
}

get_diseasestate_p_dom<-function(burdenresult,P_DOM,CASE_TOTAL_AC){
  disease<-c()
  for(i in 1:length(burdenresult$entrezgene_id)){
    if((burdenresult[i,"P_DOM"]<P_DOM)|(burdenresult[i,"CASE_TOTAL_AC"]>CASE_TOTAL_AC)){
      disease<-c(disease,1)
    }else{disease<-c(disease,0)}
  }
  return(disease)
}

get_module_pathway<-function(module_list,geneinfo){
  enrichBP_module<-data.frame()
  for(i in 1:length(module_list)){
    module_df<-geneinfo[which(geneinfo$ModuleColor==module_list[i]),]
    module_gene<-module_df$EntrezID
    enrichBP_tmp<-enrichGO(module_gene,OrgDb="org.Hs.eg.db",keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.1,pAdjustMethod = "BH",qvalueCutoff = 0.05,minGSSize = 10,maxGSSize = 500,readable = FALSE,pool = FALSE)
    enrichBP_tmp<-as.data.frame(enrichBP_tmp)
    enrichBP_tmp$modulecolor<-rep(module_list[i],length(enrichBP_tmp[,1]))
    enrichBP_module<-rbind(enrichBP_module,enrichBP_tmp)
  }
  return(enrichBP_module)
}

get_module_pathway_gprofiler<-function(module_list,geneinfo){
  enrichBP_module<-data.frame()
  for(i in 1:length(module_list)){
    module_df<-geneinfo[which(geneinfo$ModuleColor==module_list[i]),]
    module_gene<-module_df$EntrezID
    enrichBP_tmp<-gost(query = module_gene,
                                 organism = "hsapiens", ordered_query = FALSE, 
                                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                 measure_underrepresentation = FALSE, evcodes = TRUE,
                                 user_threshold = 0.1, correction_method = "gSCS", 
                                 domain_scope = "annotated", custom_bg = NULL,
                                 numeric_ns = "", sources = "GO:BP")
    enrichBP_tmp_output<-as.data.frame(enrichBP_tmp[["result"]])
    if(length(enrichBP_tmp[["result"]])!=0){
      enrichBP_tmp_output<-enrichBP_tmp_output[,which(names(enrichBP_tmp_output)!="parents")]
      enrichBP_tmp_output$modulecolor<-rep(module_list[i],length(enrichBP_tmp_output[,1]))
      enrichBP_module<-rbind(enrichBP_module,enrichBP_tmp_output)
    }
  }
  return(enrichBP_module)
}

get_heatmatrix<-function(module_enrich,all_module_enrich_names){
  class_list<-c()
  group_class<-unique(module_enrich$group)
  for(i in 1:length(group_class)){
    module_enrich_group<-module_enrich[which(module_enrich$group==group_class[i]),]
    for(j in 1:length(all_module_enrich_names$Var1)){
      if(all_module_enrich_names[j,"Var1"] %in% module_enrich_group[,"modulecolor"]){
        class_list<-c(class_list,module_enrich_group[which(module_enrich_group$modulecolor==all_module_enrich_names[j,"Var1"]),"p_adjust"])
      }else{class_list<-c(class_list,"NA")}
    }
  }
  return(class_list)
}

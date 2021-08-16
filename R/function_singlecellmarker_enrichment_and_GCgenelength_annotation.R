get_cellmarker_identity<-function(burdenresult,celltype_markers){
  #adding the title using the celltype marker list
  celltype_names<-unique(celltype_markers$cell.type)
  celltype_names<-gsub(" ","_",celltype_names)
  celltype_df<-burdenresult$GENE
  #build 0 matrix using the matrix
  for(i in 1:length(celltype_names)){
    assign(celltype_names[i],rep(0,length(burdenresult[,1])))
    celltype_df<-cbind(celltype_df,get(celltype_names[i]))
  }
  colnames(celltype_df)<-c("GENE",celltype_names)
  burdenresult<-merge(burdenresult, celltype_df, by = "GENE")
  return(burdenresult)
}

get_cellmarker_annotation<-function(burdenresult,celltype_list,celltype_markers){
  for(i in 1:length(celltype_list)){
    cellmarker_df<-celltype_markers[which(celltype_markers$cell.type==celltype_list[i]),]
    burdenresult[,16+i]<-unlist(lapply(burdenresult$entrez_id,function(x){if(x %in% cellmarker_df$entrez_id){1}else{0}}))
  }
  return(burdenresult)
}

get_genelengthGCcontent<-function(burdenresult){
  # dependent:
  #library("EDASeq")
  #library("biomaRt")
  #mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  burdenresult_list_all<-burdenresult$entrez_id
  ensembl_list_all<-getBM(values = burdenresult_list_all, filters = "entrezgene_id", mart= mart, attributes = c("entrezgene_id","ensembl_gene_id","percentage_gene_gc_content"))
  GeneLength<-getBM(values = ensembl_list_all$ensembl_gene_id, filters = "ensembl_gene_id", mart= mart, attributes = c("ensembl_gene_id","transcript_length"))
  GeneLength<-by(GeneLength,GeneLength$ensembl_gene_id,function(x){x[which.max(x$transcript_length),"transcript_length"]})
  GeneLength_df<-cbind(ensembl_gene_id=names(GeneLength),transcript_length=as.integer(GeneLength))
  GeneLength_df<-as.data.frame(GeneLength_df)
  burdengene_character<-merge(ensembl_list_all,GeneLength_df,by="ensembl_gene_id")
  burdengene_character<-burdengene_character[!duplicated(burdengene_character$ensembl_gene_id),]
  return(burdengene_character)
}

get_celltypeEnrichment<-function(burdenresult){
  all_fit_disease<-data.frame()
  burdenresult_colnames<-colnames(burdenresult)
  for(i in 20:length(burdenresult[1,])){
    burdenresult.fit<-glm(burdenresult[,i]~disease+percentage_gene_gc_content+transcript_length,data = burdenresult)
    tmp<-summary(burdenresult.fit)
    tmp<-as.data.frame(tmp$coefficients)
    tmp<-tmp["disease",]
    rownames(tmp)<-burdenresult_colnames[i]
    all_fit_disease<-rbind(all_fit_disease,tmp)
  }
  all_fit_disease$cell_type<-rownames(all_fit_disease)
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

get_heatmatrix<-function(celltype_enrich,all_celltype_enrich_names){
  class_list<-c()
  group_class<-unique(celltype_enrich$group)
  for(i in 1:length(group_class)){
    celltype_enrich_group<-celltype_enrich[which(celltype_enrich$group==group_class[i]),]
    for(j in 1:length(all_celltype_enrich_names$Var1)){
      if(all_celltype_enrich_names[j,"Var1"] %in% celltype_enrich_group[,"cell_type"]){
        class_list<-c(class_list,celltype_enrich_group[which(celltype_enrich_group$cell_type==all_celltype_enrich_names[j,"Var1"]),"p_adjust"])
      }else{class_list<-c(class_list,"NA")}
    }
  }
  return(class_list)
}

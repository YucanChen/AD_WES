#-------------------------------------------------------------------------------
# PanglaoDB database: cell type marker enrichment testing
#-------------------------------------------------------------------------------
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
  ensembl_list_all<-getBM(values = burdenresult_list_all, filters = "entrezgene_id", mart= mart, attributes = c("entrezgene_id","ensembl_gene_id","transcript_length"))
  ensembl_list_all_by<-by(ensembl_list_all,ensembl_list_all$entrezgene_id,function(x){x[which.max(x$transcript_length),c("ensembl_gene_id","transcript_length")]})
  tmp<-unlist(ensembl_list_all_by)
  tmp<-as.data.frame(matrix(tmp,nrow=length(names(ensembl_list_all_by)),ncol = 2,byrow=TRUE))
  colnames(tmp)<-c("ensembl_gene_id","transcript_length")
  ensembl_list_all<-cbind(entrezgene_id=names(ensembl_list_all_by),tmp)
  GCcontent<-getBM(values = ensembl_list_all$ensembl_gene_id, filters = "ensembl_gene_id", mart= mart, attributes = c("ensembl_gene_id","percentage_gene_gc_content"))
  burdengene_character<-merge(ensembl_list_all,GCcontent,by=c("ensembl_gene_id"))
  burdengene_character$transcript_length<-as.integer(burdengene_character$transcript_length)
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

#-------------------------------------------------------------------------------
# CellMarker database: cell type marker enrichment testing
#-------------------------------------------------------------------------------
get_CellMarker_celltypefilter<-function(celltype_marker,filtertype,filter_value){
  # <dependencies>: dplyr,magrittr,tidyverse
  if("speciesType" %in% filtertype){
    celltype_marker<-filter(celltype_marker,speciesType==filter_value[filtertype %in% "speciesType"])
  }
  if("tissueType" %in% filtertype){
    celltype_marker<-filter(celltype_marker,tissueType==filter_value[filtertype %in% "tissueType"])
  }
  if("cancerType" %in% filtertype){
    celltype_marker<-filter(celltype_marker,cancerType==filter_value[filtertype %in% "cancerType"])
  }
  celltype_marker$tissue_cellname <- paste(celltype_marker$tissueType,celltype_marker$cellName,sep = '__')
  celltype_marker_new <- celltype_marker %>%
    group_by(speciesType,tissueType,cancerType,cellType,cellName,CellOntologyID,tissue_cellname) %>%
    summarise(geneSymbol=paste(gsub("\\[|\\]","",geneSymbol),collapse=", "),geneID=paste(gsub("\\[|\\]","",geneID),collapse=", ")) %>%
    apply(1, function(x){
      x<-as.data.frame(t(x))
      gene_df<-x %>% dplyr::select(geneSymbol,geneID)
      tmp.geneID<-gene_df$geneID %>% strsplit(", ") %>% unlist %>% as.data.frame %>% dplyr::rename(c(geneID=.)) %>% filter(geneID!="NA")
      x %<>% dplyr::select(-c("geneSymbol","geneID"))
      gene_df<-cbind(x[rep(1:nrow(x),nrow(tmp.geneID)),],tmp.geneID)
      return(gene_df)
    })
  celltype_marker_new <- Reduce(function(x,y) merge(x,y,all=T),celltype_marker_new)
  return(celltype_marker_new)
}

get_CellMarker_annotation_prepare<-function(df.geneID_tissuecellname,background_genes){
  # use the TRAPD output
  annotate_df<-background_genes %>%
               dplyr::select(entrezgene_id,ensembl_gene_id,transcript_length,percentage_gene_gc_content) %>%
               group_by(entrezgene_id,ensembl_gene_id,transcript_length,percentage_gene_gc_content) %>%
               summarise()
  celltype_list<-unique(df.geneID_tissuecellname$tissue_cellname)
  col_num<-ncol(annotate_df)
  for(i in 1:length(celltype_list)){
    celltype_list_genes<-df.geneID_tissuecellname %>% filter(tissue_cellname==celltype_list[i])
    annotate_df[,col_num+i]<-unlist(lapply(annotate_df$entrezgene_id,function(x){if(x %in% celltype_list_genes$geneID){1}else{0}}))
  }
  colnames(annotate_df)<-c(colnames(annotate_df)[1:col_num],celltype_list)
  annotate_df %<>% data.frame
  return(annotate_df)
}

get_CellMarker_celltypeEnrichment<-function(annotate_df,disease_genes){
  # input disease_genes as character
  annotate_df$disease<-unlist(lapply(annotate_df$entrezgene_id,function(x){if(x %in% disease_genes){1}else{0}}))
  all_fit_disease<-data.frame()
  annotate_df.colnames<-colnames(annotate_df)
  for(i in (4+1):(ncol(annotate_df)-1)){
    annotate_df.fit<-glm(as.numeric(unlist(annotate_df[,i]))~disease+percentage_gene_gc_content+transcript_length,data = annotate_df)
    tmp<-summary(annotate_df.fit)$coefficients %>% as.data.frame
    tmp<-tmp["disease",]
    rownames(tmp)<-annotate_df.colnames[i]
    all_fit_disease<-rbind(all_fit_disease,tmp)
  }
  all_fit_disease$cell_type<-rownames(all_fit_disease)
  all_fit_disease$p_adjust<-p.adjust(all_fit_disease$`Pr(>|t|)`,method = "fdr")
  all_fit_disease<-all_fit_disease[order(all_fit_disease$p_adjust),]
  return(all_fit_disease)
}

get_CellMarker_enrichfilter<-function(df.CellMarker_enrich_res,p_adjust_value){
  df.CellMarker_enrich_res %<>% filter(p_adjust<p_adjust_value)
  if(nrow(df.CellMarker_enrich_res)>0){
    return(df.CellMarker_enrich_res)
  }
  else{
    print("Not exist!")
  }
}

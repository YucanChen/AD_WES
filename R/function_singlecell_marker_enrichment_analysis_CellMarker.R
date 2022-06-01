# annotate the background genes with gene length and GC content info

#' @import biomaRt
#' mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

get_genelengthGCcontent<-function(burdenresult){
  burdenresult_list_all<-burdenresult$entrezgene_id
  ensembl_list_all<-getBM(values = burdenresult_list_all, filters = "entrezgene_id", mart= mart, attributes = c("entrezgene_id","ensembl_gene_id","transcript_length"))
  # annotate transcript_length
  ensembl_list_all_by<-by(ensembl_list_all,ensembl_list_all$entrezgene_id,function(x){x[which.max(x$transcript_length),c("ensembl_gene_id","transcript_length")]})
  tmp<-unlist(ensembl_list_all_by) %>% matrix(nrow=length(names(ensembl_list_all_by)),ncol = 2,byrow=TRUE) %>% as.data.frame
  colnames(tmp)<-c("ensembl_gene_id","transcript_length")
  ensembl_list_all<-cbind(entrezgene_id=names(ensembl_list_all_by),tmp)
  # annotate percentage_gene_gc_content
  GCcontent<-getBM(values = ensembl_list_all$ensembl_gene_id, filters = "ensembl_gene_id", mart= mart, attributes = c("ensembl_gene_id","percentage_gene_gc_content"))
  burdengene_character<-merge(ensembl_list_all,GCcontent,by=c("ensembl_gene_id"))
  burdengene_character$transcript_length %<>% as.integer
  burdengene_character %<>% filter(!duplicated(ensembl_gene_id))
  # merge table for output
  burdengene_character %<>% merge(burdenresult,by=c("entrezgene_id")) %>%
    group_by(entrezgene_id,ensembl_gene_id,GENE,transcript_length,percentage_gene_gc_content,P_DOM_adjustFDR,disease) %>%
    summarise() %>% as.data.frame
  return(burdengene_character)
}

#-------------------------------------------------------------------------------
# CellMarker: cell type marker enrichment testing
#-------------------------------------------------------------------------------

#' @import magrittr,tidyverse

# filter the cell marker genes

get_CellMarker_celltypefilter<-function(celltype_marker,filtertype,filter_value){
  if("speciesType" %in% filtertype){
    celltype_marker<-filter(celltype_marker,speciesType==filter_value[filtertype %in% "speciesType"])
  }
  if("tissueType" %in% filtertype){
    celltype_marker<-filter(celltype_marker,tissueType==filter_value[filtertype %in% "tissueType"])
  }
  if("cancerType" %in% filtertype){
    celltype_marker<-filter(celltype_marker,cancerType==filter_value[filtertype %in% "cancerType"])
  }
  # from cellMarker delete cellMarker, proteinName, proteinID, markerResource, PMID, Company
  celltype_marker$tissue_cellname <- paste(celltype_marker$tissueType,celltype_marker$cellName,sep = '__')
  # combine table according to speciesType, tissueType, cancerType, cellType, cellName
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

# annotate the background gene tables with cell markers for enrichment test

get_CellMarker_annotation_prepare<-function(df.geneID_tissuecellname,background_genes){
  # use the preprocessed background gene tables from the TRAPD outputs
  celltype_list<-unique(df.geneID_tissuecellname$tissue_cellname)
  col_num<-ncol(background_genes)
  for(i in 1:length(celltype_list)){
    celltype_list_genes<-df.geneID_tissuecellname %>% filter(tissue_cellname==celltype_list[i])
    background_genes[,col_num+i]<-unlist(lapply(background_genes$entrezgene_id,function(x){if(x %in% celltype_list_genes$geneID){1}else{0}}))
  }
  colnames(background_genes)<-c(colnames(background_genes)[1:col_num],celltype_list)
  background_genes %<>% data.frame
  return(background_genes)
}

# conduct enrichment test to the candidate gene lists

get_CellMarker_celltypeEnrichment<-function(annotate_df,disease_genes){
  # input disease_genes as character
  annotate_df$disease<-unlist(lapply(annotate_df$entrezgene_id,function(x){if(x %in% disease_genes){1}else{0}}))
  all_fit_disease<-data.frame()
  annotate_df.colnames<-colnames(annotate_df)
  for(i in (7+1):ncol(annotate_df)){
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

# filter the significant results according to the set threshold

get_CellMarker_enrichfilter<-function(df.CellMarker_enrich_res,p_adjust_value){
  df.CellMarker_enrich_res %<>% filter(p_adjust<p_adjust_value)
  if(nrow(df.CellMarker_enrich_res)>0){
    print("exist!")
    return(df.CellMarker_enrich_res)
  }
  else{
    print("Not exist!")
  }
}

#-------------------------------------------------------------------------------
# Plot drawing
#-------------------------------------------------------------------------------
#' @import lattice,RColorBrewer

get_heatmatrix<-function(celltype_enrich,ls.sig_CellMarker,testing_set_name){
  # prepare padj value table & estimate value table
  class_list.padj<-c()
  class_list.estimate<-c()
  group_class<-unique(celltype_enrich$group)
  for(i in 1:length(group_class)){
    celltype_enrich_group<-celltype_enrich %>% filter(group==group_class[i])
    for(j in 1:length(ls.sig_CellMarker)){
      if(ls.sig_CellMarker[j] %in% celltype_enrich_group[,"cell_type"]){
        # padj value
        class_list.padj<-class_list.padj %<>% c(celltype_enrich_group[which(celltype_enrich_group$cell_type==ls.sig_CellMarker[j]),"p_adjust"])
        # Estimate value
        est_value<-celltype_enrich_group[which(celltype_enrich_group$cell_type==ls.sig_CellMarker[j]),"Estimate"]
        class_list.estimate<-class_list.estimate %<>% 
          c(if(est_value>0){1}else if(est_value<0){-1}else{0})
      }
      else{
        class_list.padj %<>% c("NA")
        class_list.estimate %<>% c("NA")
        }
    }
  }
  # prepare plotting table
  y.padj<-class_list.padj %>% matrix(., ncol = length(ls.sig_CellMarker), byrow = TRUE,dimnames=list(group_class,ls.sig_CellMarker))
  y.padj %<>% log10 %>% multiply_by(-1)
  y.estimate<-class_list.estimate %>% matrix(., ncol = length(ls.sig_CellMarker), byrow = TRUE,dimnames=list(group_class,ls.sig_CellMarker))
  y.plt<-y.padj*y.estimate
  
  # mark significant status in dataframe: y.plt
  plt.sig_statu<-y.padj
  plt.sig_statu[y.padj <= -log10(0.1)] <- ""
  plt.sig_statu[y.padj > -log10(0.001)] <- "***"
  plt.sig_statu[(y.padj > -log10(0.05))&(y.padj <= -log10(0.001))] <- "**"
  plt.sig_statu[(y.padj > -log10(0.1))&(y.padj <= -log10(0.05))] <- "*"
  
  #plt.sig_statu <- y.padj > -log10(0.1)
  #plt.sig_statu[plt.sig_statu==FALSE] <- ""
  #plt.sig_statu[plt.sig_statu==TRUE] <- "*"
  
  # visualize
  coul <- colorRampPalette(brewer.pal(8, "RdPu"))(25)
  dat <- data.frame(expand.grid(x = colnames(y.plt),y = rownames(y.plt)),value=as.numeric(t(y.plt)),sig_statu=as.character(t(plt.sig_statu)))
  Obj <- 
    levelplot(value ~ x+y,data = dat,col.regions=coul,xlab = list("cell types",fontface='bold'),ylab = list("groups",fontface='bold'),scales=list(x=list(rot=45)),main = testing_set_name) + 
    xyplot(y ~ x, data = dat,
           panel = function(y, x, ...){
             ltext(x = x, y = y, labels = dat$sig_statu, cex = 1.1, font = 1,col = '#333333')
           })
  print(Obj)
}

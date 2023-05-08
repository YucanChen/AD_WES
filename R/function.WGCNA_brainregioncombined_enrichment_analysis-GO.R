#' @import operators,tidyverse,magrittr,WGCNA,stringr

# filter the significantly correlated modules

get_sig_trait_Module<-function(moduleTraitPvalue,moduleTraitCor,set_num,trait,cutoff_pval,largerorsmaller){
  df.pval<-moduleTraitPvalue[[set_num]] %>% as.data.frame %>% dplyr::select(all_of(trait))
  df.cor<-moduleTraitCor[[set_num]] %>% as.data.frame %>% dplyr::select(all_of(trait))
  colnames(df.pval) %<>% paste(.,'.pval',sep = '')
  colnames(df.cor) %<>% paste(.,'.cor',sep = '')
  df<-cbind(df.pval,df.cor)
  df$MEColorNames<-labels2colors(as.numeric(gsub("ME","",rownames(df))))
  for(i in 1:length(trait)){
    if(largerorsmaller[i]==">"){
      df %<>% dplyr::filter(get(paste(trait[i],'.pval',sep = ''))>cutoff_pval)
    }else if(largerorsmaller[i]==">="){
      df %<>% dplyr::filter(get(paste(trait[i],'.pval',sep = ''))>=cutoff_pval)
    }else if(largerorsmaller[i]=="<"){
      df %<>% dplyr::filter(get(paste(trait[i],'.pval',sep = ''))<cutoff_pval)
    }else if(largerorsmaller[i]=="<="){
      df %<>% dplyr::filter(get(paste(trait[i],'.pval',sep = ''))<=cutoff_pval)
    }
  }
  if(nrow(df)!=0){
    return(df)
  }else{
    print("not exist!")
  }
}


#-------------------------------------------------------------------------------
# WGCNA module enrichment testing
#-------------------------------------------------------------------------------

# annotate the background gene tables with module info for enrichment test

get_WGCNAmodule_annotation<-function(df.background,geneinfo){
  # adding the title using the module gene list
  # delete grey module gene info from the gene set
  module_names<-unique(geneinfo$ModuleColor)
  module_names<-module_names[module_names %!in% "grey"]
  module_df<-df.background$entrezgene_id
  # prepare matrix for testing
  for(i in 1:length(module_names)){
    module_genes<-geneinfo %>% dplyr::filter(ModuleColor==module_names[i]) %>% dplyr::select("EntrezID")
    module_col<-df.background$entrezgene_id %>% lapply(.,function(x){if(x %in% module_genes$EntrezID){1}else{0}}) %>% unlist
    module_df<-cbind(module_df,module_col)
  }
  module_df %<>% as.data.frame %>% lapply(.,as.numeric) %>% as.data.frame
  colnames(module_df)<-c("entrezgene_id",module_names)
  df.background<-merge(df.background, module_df, by = "entrezgene_id")
  return(df.background)
}

# conduct enrichment test to the candidate gene lists

get_WGCNAmodule_Enrichment<-function(annotate_df,disease_genes){
  # mark genes with disease state (0/1)
  annotate_df %<>% tibble::add_column(.,disease=annotate_df$entrezgene_id %>% lapply(.,function(x){if(x %in% disease_genes){1}else{0}}) %>% unlist,.after = 'percentage_gene_gc_content')
  all_fit_disease<-data.frame()
  annotate_df.colnames<-colnames(annotate_df)
  # enrichment test
  for(i in (6+1):ncol(annotate_df)){
    annotate_df.fit<-glm(as.numeric(unlist(annotate_df[,i]))~disease+percentage_gene_gc_content+transcript_length,data = annotate_df)
    tmp<-summary(annotate_df.fit)$coefficients %>% as.data.frame
    tmp<-tmp["disease",]
    rownames(tmp)<-annotate_df.colnames[i]
    all_fit_disease<-rbind(all_fit_disease,tmp)
  }
  all_fit_disease$ModuleColor<-rownames(all_fit_disease)
  all_fit_disease$p_adjust<-p.adjust(all_fit_disease$`Pr(>|t|)`,method = "fdr")
  all_fit_disease<-all_fit_disease[order(all_fit_disease$p_adjust),]
  return(all_fit_disease)
}

# filter the significant results according to the set threshold

get_WGCNAmodule_enrichfilter<-function(df.WGCNAmodule_enrich_res,p_adjust_value){
  df.WGCNAmodule_enrich_res %<>% filter(p_adjust<p_adjust_value)
  if(nrow(df.WGCNAmodule_enrich_res)>0){
    print("exist!")
    return(df.WGCNAmodule_enrich_res)
  }else{
    print("Not exist!")
  }
}


#-------------------------------------------------------------------------------
# module enrichment Plot drawing
#-------------------------------------------------------------------------------
#' @import lattice,RColorBrewer

get_heatmatrix<-function(WGCNAmodule_enrich,ls.sig_WGCNAmodule,testing_set_name){
  # prepare padj value table & estimate value table
  class_list.padj<-c()
  class_list.estimate<-c()
  group_class<-unique(WGCNAmodule_enrich$group)
  for(i in 1:length(group_class)){
    WGCNAmodule_enrich_group<-WGCNAmodule_enrich %>% dplyr::filter(group==group_class[i])
    for(j in 1:length(ls.sig_WGCNAmodule)){
      if(ls.sig_WGCNAmodule[j] %in% WGCNAmodule_enrich_group[,"ModuleColor"]){
        # padj value
        class_list.padj<-class_list.padj %<>% c(WGCNAmodule_enrich_group[which(WGCNAmodule_enrich_group$ModuleColor==ls.sig_WGCNAmodule[j]),"p_adjust"])
        # Estimate value
        est_value<-WGCNAmodule_enrich_group[which(WGCNAmodule_enrich_group$ModuleColor==ls.sig_WGCNAmodule[j]),"Estimate"]
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
  y.padj<-class_list.padj %>% matrix(., ncol = length(ls.sig_WGCNAmodule), byrow = TRUE,dimnames=list(group_class,ls.sig_WGCNAmodule))
  y.padj %<>% log10 %>% multiply_by(-1)
  y.estimate<-class_list.estimate %>% matrix(., ncol = length(ls.sig_WGCNAmodule), byrow = TRUE,dimnames=list(group_class,ls.sig_WGCNAmodule))
  y.plt<-y.padj*y.estimate
  # mark significant status in dataframe: y.plt
  plt.sig_statu<-y.padj
  plt.sig_statu[y.padj <= -log10(0.1)] <- ""
  plt.sig_statu[y.padj > -log10(0.001)] <- "***"
  plt.sig_statu[(y.padj > -log10(0.05))&(y.padj <= -log10(0.001))] <- "**"
  plt.sig_statu[(y.padj > -log10(0.1))&(y.padj <= -log10(0.05))] <- "*"
  # visualize
  coul <- colorRampPalette(brewer.pal(8, "BuPu"))(25)
  dat <- data.frame(expand.grid(x = colnames(y.plt),y = rownames(y.plt)),value=as.numeric(t(y.plt)),sig_statu=as.character(t(plt.sig_statu)))
  Obj <- 
    levelplot(value ~ x+y,data = dat,col.regions=coul,xlab = list("ModuleColor",fontface='bold'),ylab = list("groups",fontface='bold'),scales=list(x=list(rot=45)),main = testing_set_name) + 
    xyplot(y ~ x, data = dat,
           panel = function(y, x, ...){
             ltext(x = x, y = y, labels = dat$sig_statu, cex = 1.1, font = 1,col = '#333333')
           })
  print(Obj)
}


#-------------------------------------------------------------------------------
# GOBP enrichment test & Plot drawing
#-------------------------------------------------------------------------------
#' @import gprofiler2,clusterProfiler,org.Hs.eg.db

# enrichment test

## gprofiler

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
      enrichBP_tmp_output$ModuleColor<-rep(module_list[i],nrow(enrichBP_tmp_output))
      enrichBP_module<-rbind(enrichBP_module,enrichBP_tmp_output)
    }
  }
  return(enrichBP_module)
}

## clusterprofiler

get_module_pathway_clusterprofiler<-function(module_list,geneinfo){
  enrichBP_module<-data.frame()
  for(i in 1:length(module_list)){
    module_df<-geneinfo[which(geneinfo$ModuleColor==module_list[i]),]
    module_gene<-module_df$EntrezID
    enrichBP_tmp<-enrichGO(module_gene,OrgDb="org.Hs.eg.db",keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.1,pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500,readable = FALSE,pool = FALSE)
    enrichBP_tmp<-as.data.frame(enrichBP_tmp)
    enrichBP_tmp$ModuleColor<-rep(module_list[i],nrow(enrichBP_tmp))
    enrichBP_module<-rbind(enrichBP_module,enrichBP_tmp)
  }
  return(enrichBP_module)
}


# enriched pathway plot drawing

get_pathway_barplot<-function(module_pathway_clusterprofiler,module.list,term_num.list){
  
  # filter modules & terms
  
  res.plt<-data.frame()
  for(i in 1:length(module.list)){
    module.terms<-module_pathway_clusterprofiler %>% dplyr::filter(ModuleColor==module.list[i])
    if(nrow(module.terms)>=term_num.list[i]){
      module.terms %<>% slice_head(n = term_num.list[i])
      res.plt %<>% rbind(.,module.terms)
    }else{
      res.plt %<>% rbind(.,module.terms)
    }
  }
  
  # barplot
  
  plt<-list()
  for(i in 1:length(module.list)){
    module.plt<-res.plt %>% dplyr::filter(ModuleColor==module.list[i])
    plt.tmp<-
      module.plt %>%
      mutate(Description = fct_reorder(Description, -log10(p.adjust))) %>%
      ggplot(aes(x=Description,y= -log10(p.adjust),fill=Count)) +
      geom_bar(stat='identity',position = "dodge",width=0.9) +
      scale_fill_gradient(expression(Count),low="#6666FF", high = "#660099") +
      coord_flip() +
      scale_x_discrete(labels=function(x) str_wrap(x, width=35)) +
      labs(x='',y='-log10(p.adjust)',title=paste('top:',nrow(module.plt),'enriched','in',module.list[i],sep = ' ')) +
      theme(
        axis.text.x=element_text(color="black",size=rel(1.1)),
        axis.text.y=element_text(color="black", size=rel(1.3)),
        axis.title.x = element_text(color="black", size=rel(1.1)),
        legend.text=element_text(color="black",size=rel(0.9)),
        legend.title = element_text(color="black",size=rel(1.0))
      ) +
      geom_hline(yintercept = -log10(0.1),color = "grey",linetype='dashed')+
      theme(panel.background = element_rect(color = 'grey60', fill = 'transparent'))
    plt[[i]]<-plt.tmp
  }
  
  return(plt)
}

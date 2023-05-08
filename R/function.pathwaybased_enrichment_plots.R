# pathway enrichment analysis

#' @import tidyverse,magrittr,gprofiler2,clusterProfiler,org.Hs.eg.db,RColorBrewer,stringr

## gprofiler

get_pathway_gprofiler<-function(gene_lists,genename_lists,sub_ontology,padj_value){
  enrich_category_table<-data.frame()
  for(i in 1:length(gene_lists)){
    list_genes<-gene_lists[i][[1]]
    enrich_category_tmp<-gost(query = list_genes,
                       organism = "hsapiens", ordered_query = FALSE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = TRUE,
                       user_threshold = padj_value, correction_method = "gSCS", 
                       domain_scope = "annotated", custom_bg = NULL,
                       numeric_ns = "", sources = paste("GO",sub_ontology,sep = ':'))
    enrich_category_tmp_output<-as.data.frame(enrich_category_tmp[["result"]])
    if(length(enrich_category_tmp[["result"]])!=0){
      enrich_category_tmp_output<-enrich_category_tmp_output[,which(names(enrich_category_tmp_output)!="parents")]
      enrich_category_tmp_output$group<-rep(genename_lists[i],nrow(enrich_category_tmp_output))
      enrich_category_table<-rbind(enrich_category_table,enrich_category_tmp_output)
    }
  }
  return(enrich_category_table)
}

## clusterprofiler

get_pathway_clusterprofiler<-function(gene_lists,genename_lists,sub_ontology,padj_value){
  enrich_category_table<-data.frame()
  for(i in 1:length(gene_lists)){
    list_genes<-gene_lists[i][[1]]
    enrich_category_tmp<-enrichGO(list_genes,OrgDb="org.Hs.eg.db",keyType = "ENTREZID",ont = sub_ontology, pvalueCutoff = padj_value, pAdjustMethod = "BH",minGSSize = 10,maxGSSize = 500,readable = FALSE,pool = FALSE)
    enrich_category_tmp<-as.data.frame(enrich_category_tmp)
    enrich_category_tmp$group<-rep(genename_lists[i],nrow(enrich_category_tmp))
    enrich_category_table<-rbind(enrich_category_table,enrich_category_tmp)
  }
  return(enrich_category_table)
}

# plot drawing: grouping bar chart

get_GOres_bubblechart<-function(GOres_combined,term_num.list,out_group.list,palette_name="BrBG",plot_title){
  
  # prepare data
  
  group.list<-out_group.list
  group.df<-data.frame()
  for(i in 1:length(group.list)){
    group.tmp<-GOres_combined %>% dplyr::filter(group==group.list[i])
    if(nrow(group.tmp)>term_num.list[i]){
      group.tmp %<>% slice_head(n = term_num.list[i])
    }
    group.df<-rbind(group.df,group.tmp)
  }
  group.df$LP <- -log10(group.df$p.adjust)
  
  # theme & color palette setting
  
  group.df <- within(group.df,{bubblesize<- sqrt((group.df$Count+100)/(100*pi))})
  group.df$Description<-str_wrap(group.df$Description,25)
  group.df$Description<-factor(group.df$Description,levels=unique(group.df$Description))
  LP.cutoff <- -log10(0.05)
  y.max <- max(group.df$LP+group.df$bubblesize)
  mytheme <- theme_minimal()+
    theme(
      panel.grid.major.y=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1,color = 'black'),
      plot.title=element_text(hjust =0.5),
      axis.line.y=element_line(linetype=1,color='grey'),
      axis.line.x=element_line(linetype=1,color='grey'),
      axis.ticks = element_line(linetype=2,color='grey'),
      panel.grid=element_line(linetype=2,color='grey'),
      legend.background = element_rect(fill="gray90", size=0,color='white'),
      legend.text=element_text(face="plain",size=11,color = 'black'),
      legend.title=element_text(face="bold",size=11,hjust = 0.05,color = 'black'),
      axis.text=element_text(face="plain",size=10)
    )
  coul <- brewer.pal(n=length(group.list), palette_name)
  
  # draw plot
  
  ggplot(data=group.df, mapping=aes(x=Description,y=LP,color=factor(group,levels=unique(group))))+
    geom_point(stat= "identity",aes(size=bubblesize),alpha=0.7,show.legend = TRUE)+ 
    guides(color=guide_legend(title="group"))+
    scale_size(range = c(1, 30),guide="none")+
    scale_color_manual(values=coul)+
    scale_y_continuous(limits=c(1,y.max))+
    labs(x='terms',y='-log10(p.adjust)',title=plot_title)+
    mytheme+
    geom_hline(yintercept = LP.cutoff,linetype='dashed')
  
}


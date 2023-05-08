deduplicate_SYMBOL <- function(count_txt){
  count_txt$SYMBOL %<>% as.factor()
  count_txt.count <- with(count_txt,
                          by(count_txt,SYMBOL,function(x){mean(x$Counts)})
                          ) %>% cbind %>% data.frame()
  count_new<-data.frame(SYMBOL=rownames(count_txt.count),count_txt.count %>% rename(.,Counts=.))
  rm(count_txt,count_txt.count)
  return(count_new)
}

mutationratio_test <- function(out_folder,PAV_count.case,PAV_count.control,Syn_count.case,Syn_count.control,merge_method,count_adjust_method,correction_method,padjust.threshold){
  
  PAV_count.case %<>% rename(.,PAV_count.case=Counts)
  PAV_count.control %<>% rename(.,PAV_count.control=Counts)
  Syn_count.case %<>% rename(.,Syn_count.case=Counts)
  Syn_count.control %<>% rename(.,Syn_count.control=Counts)
  
  # fisher test function definition
  
  func.fisher<-function(count_merged.PAV.adj,correction_method){
    df <- count_merged.PAV.adj %>% dplyr::select(c('PAV_count.case.adj','PAV_count.control.adj','PAV_remain.case.adj','PAV_remain.control.adj'))
    tmp<-apply(df,MARGIN = 1,function(x){
      test<-matrix(x,ncol = 2,byrow = TRUE,dimnames = list(c("gene", "remaining"),c("case", "control")))
      test.fisher<-fisher.test(test, alternative = "greater")
      c(test.fisher$p.value,test.fisher$estimate,test.fisher$conf.int[1])
    }) %>% Matrix::t() %>% as.data.frame()
    colnames(tmp)<-c('p.value','odds.ratio','conf.int.lower')
    tmp$p.adj<-p.adjust(p = tmp$p.value,method = correction_method)
    res.fisher_test<-cbind(count_merged.PAV.adj,tmp) %>% dplyr::arrange(desc(conf.int.lower)) %>% dplyr::arrange(p.adj)
    return(res.fisher_test)
  }
  
  # test all genes between case & control groups
  
  if('all genes' %in% merge_method){
    
    if('genes' %in% count_adjust_method){
      if((nrow(PAV_count.case)>0)&(nrow(PAV_count.control)>0)&(nrow(Syn_count.case)>0)&(nrow(Syn_count.control)>0)){
        
        # prepare table for testing dat[is.na(dat)] <- 0
        
        count.merged.genes.case<-merge(PAV_count.case,Syn_count.case,by = c('SYMBOL'),sort = T,all.y = T)
        count.merged.genes.case[is.na(count.merged.genes.case)] <- 0
        count.merged.genes.control<-merge(PAV_count.control,Syn_count.control,by = c('SYMBOL'),sort = T,all.y = T)
        count.merged.genes.control[is.na(count.merged.genes.control)] <- 0
        count.merged.genes.case %<>% tibble::add_column(.,PAV_count.case.adj = mapply(function(x,y){x/y},x = count.merged.genes.case$PAV_count.case,y = count.merged.genes.case$Syn_count.case))
        count.merged.genes.control %<>% tibble::add_column(.,PAV_count.control.adj = mapply(function(x,y){x/y},x = count.merged.genes.control$PAV_count.control,y = count.merged.genes.control$Syn_count.control))
        count.merged.genes.PAV.adj<-merge(count.merged.genes.case,count.merged.genes.control,by = c('SYMBOL'),sort = T,all = T)
        count.merged.genes.PAV.adj[is.na(count.merged.genes.PAV.adj)] <- 0
        rm(count.merged.genes.case,count.merged.genes.control)
        
        tmp<-c(count.merged.genes.PAV.adj$PAV_count.case.adj,count.merged.genes.PAV.adj$PAV_count.control.adj)
        tmp.zero<-which(tmp==0); tmp<-tmp[-tmp.zero]
        adj.multiple<-tmp %>% log10 %>% min %>% floor %>% abs
        count.merged.genes.PAV.adj$PAV_count.case.adj %<>% multiply_by(10^adj.multiple)
        count.merged.genes.PAV.adj$PAV_count.control.adj %<>% multiply_by(10^adj.multiple)
        sum.PAV_count.case.adj<-count.merged.genes.PAV.adj$PAV_count.case.adj %>% sum
        sum.PAV_count.control.adj<-count.merged.genes.PAV.adj$PAV_count.control.adj %>% sum
        count.merged.genes.PAV.adj %<>% tibble::add_column(.,PAV_remain.case.adj = sum.PAV_count.case.adj - count.merged.genes.PAV.adj$PAV_count.case.adj) %>% tibble::add_column(.,PAV_remain.control.adj = sum.PAV_count.control.adj - count.merged.genes.PAV.adj$PAV_count.control.adj)
        
        # fisher test
        
        res.fisher_test.genes <- func.fisher(count_merged.PAV.adj = count.merged.genes.PAV.adj,correction_method = correction_method)
        ls.res<-list()
        ls.res$all_res<-res.fisher_test.genes
        ls.res$sig_res<-res.fisher_test.genes %>% dplyr::filter(p.adj < padjust.threshold)
        write.xlsx(ls.res,paste(out_folder,'/res.fisher_test.all_gene.genes.',correction_method,'.',padjust.threshold,'.xlsx',sep = ''))
        
      }else{
        print('not comparable!')
      }
    }
    
    if('sum' %in% count_adjust_method){
      if((nrow(PAV_count.case)>0)&(nrow(PAV_count.control)>0)&(nrow(Syn_count.case)>0)&(nrow(Syn_count.control)>0)){
        
        # prepare table for testing
        
        count_merged.sum.PAV.adj<-merge(PAV_count.case,PAV_count.control,by = c('SYMBOL'),sort = T,all = T)
        count_merged.sum.PAV.adj<-merge(count_merged.sum.PAV.adj,Syn_count.case,by = c('SYMBOL'),sort = T,all.x = T)
        count_merged.sum.PAV.adj<-merge(count_merged.sum.PAV.adj,Syn_count.control,by = c('SYMBOL'),sort = T,all.x = T)
        count_merged.sum.PAV.adj[is.na(count_merged.sum.PAV.adj)] <- 0
        sum.count.PAV.case<-count_merged.sum.PAV.adj$PAV_count.case %>% sum
        sum.count.PAV.control<-count_merged.sum.PAV.adj$PAV_count.control %>% sum
        sum.count.syn.case<-count_merged.sum.PAV.adj$Syn_count.case %>% sum
        sum.count.syn.control<-count_merged.sum.PAV.adj$Syn_count.control %>% sum
        
        count_merged.sum.PAV.adj %<>% tibble::add_column(.,PAV_count.case.adj = count_merged.sum.PAV.adj$PAV_count.case %>% divide_by(sum.count.syn.case)) %>% tibble::add_column(.,PAV_count.control.adj = count_merged.sum.PAV.adj$PAV_count.control %>% divide_by(sum.count.syn.control))
        tmp<-c(count_merged.sum.PAV.adj$PAV_count.case.adj,count_merged.sum.PAV.adj$PAV_count.control.adj)
        tmp.zero<-which(tmp==0); tmp<-tmp[-tmp.zero]
        adj.multiple<-tmp %>% log10 %>% min %>% floor %>% abs
        count_merged.sum.PAV.adj$PAV_count.case.adj %<>% multiply_by(10^adj.multiple)
        count_merged.sum.PAV.adj$PAV_count.control.adj %<>% multiply_by(10^adj.multiple)
        count_merged.sum.PAV.adj %<>% tibble::add_column(.,PAV_remain.case.adj = (sum.count.PAV.case/sum.count.syn.case)*(10^adj.multiple) - count_merged.sum.PAV.adj$PAV_count.case.adj) %>% tibble::add_column(.,PAV_remain.control.adj = (sum.count.PAV.control/sum.count.syn.control)*(10^adj.multiple) - count_merged.sum.PAV.adj$PAV_count.control.adj)
        
        #sum.PAV_count.case.adj<-PAV_count.case$PAV_count.case %>% sum
        #sum.PAV_count.control.adj<-PAV_count.control$PAV_count.control %>% sum
        #count_merged.sum.PAV.adj %<>% tibble::add_column(.,PAV_remain.case.adj = sum.PAV_count.case.adj - count_merged.sum.PAV.adj$PAV_count.case.adj) %>% tibble::add_column(.,PAV_remain.control.adj = sum.PAV_count.control.adj - count_merged.sum.PAV.adj$PAV_count.control.adj)
        
        # fisher test
        
        res.fisher_test.sum <- func.fisher(count_merged.PAV.adj = count_merged.sum.PAV.adj,correction_method = correction_method)
        ls.res<-list()
        ls.res$all_res<-res.fisher_test.sum
        ls.res$sig_res<-res.fisher_test.sum %>% dplyr::filter(p.adj < padjust.threshold)
        write.xlsx(ls.res,paste(out_folder,'/res.fisher_test.all_gene.sum.',correction_method,'.',padjust.threshold,'.xlsx',sep = ''))
        
      }else{
        print('not comparable!')
      }
    } 
    
  }
  
  # test co-occurred genes between case & control groups
  
  if('co-occured genes' %in% merge_method){
    
    if('genes' %in% count_adjust_method){
      if((nrow(PAV_count.case)>0)&(nrow(PAV_count.control)>0)&(nrow(Syn_count.case)>0)&(nrow(Syn_count.control)>0)){
        
        # prepare table for testing
        
        count.merged.genes.case<-merge(PAV_count.case,Syn_count.case,by = c('SYMBOL'),sort = T,all.y = T)
        count.merged.genes.case[is.na(count.merged.genes.case)] <- 0
        count.merged.genes.control<-merge(PAV_count.control,Syn_count.control,by = c('SYMBOL'),sort = T,all.y = T)
        count.merged.genes.control[is.na(count.merged.genes.control)] <- 0
        count.merged.genes.case %<>% tibble::add_column(.,PAV_count.case.adj = mapply(function(x,y){x/y},x = count.merged.genes.case$PAV_count.case,y = count.merged.genes.case$Syn_count.case))
        count.merged.genes.control %<>% tibble::add_column(.,PAV_count.control.adj = mapply(function(x,y){x/y},x = count.merged.genes.control$PAV_count.control,y = count.merged.genes.control$Syn_count.control))
        count.merged.genes.PAV.adj<-merge(count.merged.genes.case,count.merged.genes.control,by = c('SYMBOL'),sort = T,all = F)
        rm(count.merged.genes.case,count.merged.genes.control)
        
        tmp<-c(count.merged.genes.PAV.adj$PAV_count.case.adj,count.merged.genes.PAV.adj$PAV_count.control.adj)
        tmp.zero<-which(tmp==0); tmp<-tmp[-tmp.zero]
        adj.multiple<-tmp %>% log10 %>% min %>% floor %>% abs
        count.merged.genes.PAV.adj$PAV_count.case.adj %<>% multiply_by(10^adj.multiple)
        count.merged.genes.PAV.adj$PAV_count.control.adj %<>% multiply_by(10^adj.multiple)
        sum.PAV_count.case.adj<-count.merged.genes.PAV.adj$PAV_count.case.adj %>% sum
        sum.PAV_count.control.adj<-count.merged.genes.PAV.adj$PAV_count.control.adj %>% sum
        count.merged.genes.PAV.adj %<>% tibble::add_column(.,PAV_remain.case.adj = sum.PAV_count.case.adj - count.merged.genes.PAV.adj$PAV_count.case.adj) %>% tibble::add_column(.,PAV_remain.control.adj = sum.PAV_count.control.adj - count.merged.genes.PAV.adj$PAV_count.control.adj)
        
        # fisher test
        
        res.fisher_test.genes <- func.fisher(count_merged.PAV.adj = count.merged.genes.PAV.adj,correction_method = correction_method)
        ls.res<-list()
        ls.res$all_res<-res.fisher_test.genes
        ls.res$sig_res<-res.fisher_test.genes %>% dplyr::filter(p.adj < padjust.threshold)
        write.xlsx(ls.res,paste(out_folder,'/res.fisher_test.co_occurred.genes.',correction_method,'.',padjust.threshold,'.xlsx',sep = ''))
        
      }else{
        print('not comparable!')
      }
    }
    
    if('sum' %in% count_adjust_method){
      if((nrow(PAV_count.case)>0)&(nrow(PAV_count.control)>0)&(nrow(Syn_count.case)>0)&(nrow(Syn_count.control)>0)){
        
        # prepare table for testing
        
        count_merged.sum.PAV.adj<-merge(PAV_count.case,PAV_count.control,by = c('SYMBOL'),sort = T,all = F)
        count_merged.sum.PAV.adj<-merge(count_merged.sum.PAV.adj,Syn_count.case,by = c('SYMBOL'),sort = T,all.x = T)
        count_merged.sum.PAV.adj<-merge(count_merged.sum.PAV.adj,Syn_count.control,by = c('SYMBOL'),sort = T,all.x = T)
        count_merged.sum.PAV.adj[is.na(count_merged.sum.PAV.adj)] <- 0
        sum.count.PAV.case<-count_merged.sum.PAV.adj$PAV_count.case %>% sum
        sum.count.PAV.control<-count_merged.sum.PAV.adj$PAV_count.control %>% sum
        sum.count.syn.case<-count_merged.sum.PAV.adj$Syn_count.case %>% sum
        sum.count.syn.control<-count_merged.sum.PAV.adj$Syn_count.control %>% sum
        
        count_merged.sum.PAV.adj %<>% tibble::add_column(.,PAV_count.case.adj = count_merged.sum.PAV.adj$PAV_count.case %>% divide_by(sum.count.syn.case)) %>% tibble::add_column(.,PAV_count.control.adj = count_merged.sum.PAV.adj$PAV_count.control %>% divide_by(sum.count.syn.control))
        adj.multiple<-c(count_merged.sum.PAV.adj$PAV_count.case.adj,count_merged.sum.PAV.adj$PAV_count.control.adj) %>% log10 %>% min %>% floor %>% abs
        count_merged.sum.PAV.adj$PAV_count.case.adj %<>% multiply_by(10^adj.multiple)
        count_merged.sum.PAV.adj$PAV_count.control.adj %<>% multiply_by(10^adj.multiple)
        count_merged.sum.PAV.adj %<>% tibble::add_column(.,PAV_remain.case.adj = (sum.count.PAV.case/sum.count.syn.case)*(10^adj.multiple) - count_merged.sum.PAV.adj$PAV_count.case.adj) %>% tibble::add_column(.,PAV_remain.control.adj = (sum.count.PAV.control/sum.count.syn.control)*(10^adj.multiple) - count_merged.sum.PAV.adj$PAV_count.control.adj)
        
        #sum.PAV_count.case.adj<-count_merged.sum.PAV.adj$PAV_count.case.adj %>% sum
        #sum.PAV_count.control.adj<-count_merged.sum.PAV.adj$PAV_count.control.adj %>% sum
        #count_merged.sum.PAV.adj %<>% tibble::add_column(.,PAV_remain.case.adj = sum.PAV_count.case.adj - count_merged.sum.PAV.adj$PAV_count.case.adj) %>% tibble::add_column(.,PAV_remain.control.adj = sum.PAV_count.control.adj - count_merged.sum.PAV.adj$PAV_count.control.adj)
        
        # fisher test
        
        res.fisher_test.sum <- func.fisher(count_merged.PAV.adj = count_merged.sum.PAV.adj,correction_method = correction_method)
        ls.res<-list()
        ls.res$all_res<-res.fisher_test.sum
        ls.res$sig_res<-res.fisher_test.sum %>% dplyr::filter(p.adj < padjust.threshold)
        write.xlsx(ls.res,paste(out_folder,'/res.fisher_test.co_occurred.sum.',correction_method,'.',padjust.threshold,'.xlsx',sep = ''))
        
      }else{
        print('not comparable!')
      }
    } 
    
  }
  
  print('analysis finished!')
  
}

# Task Description: Fisher’s test on mutation counts adjusted by synonymous mutation.

# 1. find the co-occurred genes between case (AD/dementia) & control groups
# 2. adjust the count values of co-occurred genes using Syn mutations:
# (1) use the total values of syn mutations from case & control respectively for adjustment;
# (2) use the gene-level values of syn mutations from case & control ......
# 3. perform gene-level fisher test for each co-occurred genes
# 4.perform fdr correction for the fisher test p-value of snp/indel/merge co-occurred genes between case & control groups

# count_adjust_method=c('sum','genes')
#out_folder<-"/AD_exome/Downstream_analyses/AD_downstream_scripts_20220627/mutationratio_validation/results/dementia"

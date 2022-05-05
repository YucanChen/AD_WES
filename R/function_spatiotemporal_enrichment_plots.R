#-------------------------------------------------------------------------------
# BrainSpan RNA-seq spatio-temporal data preprocessing
#-------------------------------------------------------------------------------

# mark the temporal development periods according to sample ages

get_temperalperiod_list <- function(Brainspan_RNAseq_col.age){
  temperal_period_list<-c()
  tmp<-c()
  for(i in 1:length(Brainspan_RNAseq_col.age$number)){
    if(Brainspan_RNAseq_col.age[i,"unit"]=="pcw"){
      if((Brainspan_RNAseq_col.age[i,"number"]>=4)&(Brainspan_RNAseq_col.age[i,"number"]<8)){tmp<-"T01_embryonic";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=8)&(Brainspan_RNAseq_col.age[i,"number"]<10)){tmp<-"T02_early_fetal1";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=10)&(Brainspan_RNAseq_col.age[i,"number"]<13)){tmp<-"T03_early_fetal2";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=13)&(Brainspan_RNAseq_col.age[i,"number"]<16)){tmp<-"T04_early_midfetal1";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=16)&(Brainspan_RNAseq_col.age[i,"number"]<19)){tmp<-"T05_early_midfetal2";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=19)&(Brainspan_RNAseq_col.age[i,"number"]<24)){tmp<-"T06_late_midfetal";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=24)&(Brainspan_RNAseq_col.age[i,"number"]<38)){tmp<-"T07_late_fetal";temperal_period_list<-c(temperal_period_list,tmp)}
      else{tmp<-c()}
    }
    
    if(Brainspan_RNAseq_col.age[i,"unit"]=="mos"){
      if((Brainspan_RNAseq_col.age[i,"number"]>=0)&(Brainspan_RNAseq_col.age[i,"number"]<6)){tmp<-"T08_neonatal_earlyinfancy";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=6)&(Brainspan_RNAseq_col.age[i,"number"]<12)){tmp<-"T09_late_infancy";temperal_period_list<-c(temperal_period_list,tmp)}
      else{tmp<-c()}
    }
    
    if(Brainspan_RNAseq_col.age[i,"unit"]=="yrs"){
      if((Brainspan_RNAseq_col.age[i,"number"]>=1)&(Brainspan_RNAseq_col.age[i,"number"]<6)){tmp<-"T10_early_childhood";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=6)&(Brainspan_RNAseq_col.age[i,"number"]<12)){tmp<-"T11_midlate_childhood";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=12)&(Brainspan_RNAseq_col.age[i,"number"]<20)){tmp<-"T12_adolescence";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=20)&(Brainspan_RNAseq_col.age[i,"number"]<40)){tmp<-"T13_young_adulthood";temperal_period_list<-c(temperal_period_list,tmp)}
      if((Brainspan_RNAseq_col.age[i,"number"]>=40)&(Brainspan_RNAseq_col.age[i,"number"]<60)){tmp<-"T14_middle_adulthood";temperal_period_list<-c(temperal_period_list,tmp)}
      else{tmp<-c()}
    }
  }
  return(temperal_period_list)
}

# obtain the spatio-temporal distribution of available samples

get_period_stucture <- function(Brainspan_RNAseq_col.age){
  period_stucture_list<-as.character(unique(Brainspan_RNAseq_col.age$structure))
  tmp<-data.frame()
  period_stucture_table <- as.data.frame(period_stucture_list)
  period_list <- unique(Brainspan_RNAseq_col.age$developmental_period)
  for(i in 1:length(period_list)){
    tmp<-as.data.frame(sort(table(Brainspan_RNAseq_col.age[which(Brainspan_RNAseq_col.age$developmental_period==period_list[i]),"structure"])))
    tmp_number_list<-c()
    for(j in 1:length(period_stucture_list)){
      tmp_number<-c()
      if(period_stucture_list[j] %in% tmp[,1]){
        tmp_number=tmp[which(tmp$Var1==period_stucture_list[j]),2]
      }else{tmp_number=0}
      tmp_number_list<-c(tmp_number_list,tmp_number)
    }
    period_stucture_table <- cbind(period_stucture_table,as.data.frame(tmp_number_list))
  }
  structure_name<-as.character(unique(Brainspan_RNAseq_col.age[(period_stucture_list %in% Brainspan_RNAseq_col.age$structure),"structure_name"]))
  period_stucture_table <- cbind(period_stucture_table,as.data.frame(structure_name))
  rownames(period_stucture_table)<-period_stucture_list
  period_stucture_table<-period_stucture_table[,-c(1)]
  colnames(period_stucture_table)<-c(period_list,"structure_name")
  return(period_stucture_table)
}

# period 3-7 data extraction for clustering (brain-region devided)

get_hierarchical_expression <- function(Brainspan_RNAseq_col.age,Brainspan_RNAseq_row,Brainspan_RNAseq){
  division_period_list<-c("T03_early_fetal2","T04_early_midfetal1","T05_early_midfetal2","T06_late_midfetal","T07_late_fetal")
  division_period_table<-data.frame(Brainspan_RNAseq_row$gene_id)
  tmp<-data.frame()
  division_period_table.col<-data.frame()
  tmp.col<-data.frame()
  for(i in 1:length(Brainspan_RNAseq_col.age$developmental_period)){
    if(Brainspan_RNAseq_col.age[i,"developmental_period"] %in% division_period_list){
      tmp<-as.data.frame(Brainspan_RNAseq[,i])
      division_period_table<-cbind(division_period_table,tmp)
      tmp.col<-as.data.frame(Brainspan_RNAseq_col.age[i,])
      division_period_table.col<-rbind(division_period_table.col,tmp.col)
    }else{tmp<-data.frame();tmp.col<-data.frame()}
  }
  division_period_table<-division_period_table[,-c(1)]
  division_period_table<-rbind(division_period_table.col$structure,division_period_table)
  division_regionselected_list<-c("Ocx","M1C-S1C","MGE","URL","CGE","DTH","LGE","PCx","TCx","CB")
  division_regionselected_table<-data.frame(c(1:length(division_period_table[,1])))
  tmp<-data.frame()
  division_regionselected_table.col<-data.frame()
  for(j in 1:length(division_period_table[1,])){
    if(!(division_period_table[1,j] %in% division_regionselected_list)){
      tmp<-as.data.frame(division_period_table[,j])
      division_regionselected_table<-cbind(division_regionselected_table,tmp)
    }else{tmp<-data.frame()}
  }
  division_regionselected_table<-division_regionselected_table[,-c(1)]
  return(division_regionselected_table)
}

# merge PCs of samples from the same brain regions via median values for hclust

get_brainregion_clustertable<-function(period3_7_table){
  region_name<-unique(as.character(unique(as.character(gsub("_.*","",colnames(period3_7_table))))))
  brainregion_clustertable<-data.frame(c("brain_region",1:length(period3_7_table[,1])))
  tmp_sum_ave<-data.frame()
  tmp_sum<-data.frame()
  tmp<-data.frame()
  period3_7_table.brain_region<-as.character(gsub("_.*","",colnames(period3_7_table)))
  for(i in 1:length(region_name)){
    tmp_sum<-as.data.frame(c(1:length(period3_7_table[,1])))
    for(j in 1:length(period3_7_table[1,])){
      if(period3_7_table.brain_region[j]==region_name[i]){
        tmp<-period3_7_table[,j]
        tmp_sum<-cbind(tmp_sum,tmp)
      }else{tmp<-data.frame()}
    }
    tmp_sum<-tmp_sum[,-c(1)]
    tmp_sum_ave<-apply(tmp_sum,1,median)
    tmp_sum_ave<-as.data.frame(c(region_name[i],tmp_sum_ave))
    colnames(brainregion_clustertable)<-c("1")
    brainregion_clustertable<-cbind(brainregion_clustertable,tmp_sum_ave)
  }
  colnames(brainregion_clustertable)<-brainregion_clustertable[1,]
  brainregion_clustertable<-brainregion_clustertable[-c(1),];brainregion_clustertable<-brainregion_clustertable[,-c(1)]
  brainregion_clustertable<-t(brainregion_clustertable)
  brainregion_clustertable.col<-colnames(brainregion_clustertable)
  brainregion_clustertable<-apply(brainregion_clustertable,1,as.numeric)
  brainregion_clustertable<-as.data.frame(t(brainregion_clustertable))
  colnames(brainregion_clustertable)<- brainregion_clustertable.col
  return(brainregion_clustertable)
}

# delete samples from specific brain regions

## delete samples from the expression matrix

deregion_Brainspan_RNAseq<-function(regionselected_list,Brainspan_RNAseq,Brainspan_RNAseq_col){
  RNAseq_regionselected_table<-data.frame(c(1:length(Brainspan_RNAseq[,1])))
  tmp<-data.frame()
  for(i in 1:length(Brainspan_RNAseq_col$structure_acronym)){
    if(!(Brainspan_RNAseq_col[i,"structure_acronym"] %in% regionselected_list)){
      tmp<-as.data.frame(Brainspan_RNAseq[,i])
      RNAseq_regionselected_table<-cbind(RNAseq_regionselected_table,tmp)
    }else{tmp<-data.frame()}
  }
  RNAseq_regionselected_table<-RNAseq_regionselected_table[,-c(1)]
  return(RNAseq_regionselected_table)
}

## delete samples from the sample column info

deregion_Brainspan_RNAseq.col<-function(regionselected_list,Brainspan_RNAseq_col){
  RNAseq_regionselected_table.col<-as.data.frame(Brainspan_RNAseq_col[1,])
  tmp<-data.frame()
  for(i in 1:length(Brainspan_RNAseq_col$structure_acronym)){
    if(!(Brainspan_RNAseq_col[i,"structure_acronym"] %in% regionselected_list)){
      tmp<-as.data.frame(Brainspan_RNAseq_col[i,])
      RNAseq_regionselected_table.col<-rbind(RNAseq_regionselected_table.col,tmp)
    }else{tmp<-data.frame()}
  }
  RNAseq_regionselected_table.col<-RNAseq_regionselected_table.col[-c(1),]
  colnames(RNAseq_regionselected_table.col)<-colnames(Brainspan_RNAseq_col)
  return(RNAseq_regionselected_table.col)
}

# mark samples with cluster info

mark_cluster_expressiontable<-function(cluster_list,Brainspan_RNAseq_col){
  n=length(cluster_list)
  tmp<-c()
  cluster_col<-c()
  for(k in 1:length(Brainspan_RNAseq_col$structure_acronym)){
    for(i in 1:n){
      for(j in 1:length(cluster_list[[i]])){
        if(Brainspan_RNAseq_col[k,"structure_acronym"]==cluster_list[[i]][j]){
          tmp<-paste("cluster",i,sep = "")
          cluster_col<-c(cluster_col,tmp)
        }else{tmp<-c()}
      }
    }
  }
  Brainspan_RNAseq_col$cluster<-cluster_col
  return(Brainspan_RNAseq_col)
}

getuniqueentrez <-function(Brainspan_RNAseq_row,Brainspan_RNAseq){
  single_entrez<-as.data.frame(table(Brainspan_RNAseq_row$entrez_id))
  single_entrez<-single_entrez[which(single_entrez$Freq==1),]
  retain_ensembl<-unlist(lapply(single_entrez$Var1, function(x){Brainspan_RNAseq_row[which(Brainspan_RNAseq_row$entrez_id==x),"ensembl_gene_id"]}))
  dup_entrez<-as.data.frame(table(Brainspan_RNAseq_row$entrez_id))
  dup_entrez<-dup_entrez[which(dup_entrez$Freq>1),]
  dup_ensembl<-Brainspan_RNAseq_row[which(Brainspan_RNAseq_row$entrez_id %in% dup_entrez$Var1),]
  #length(as.data.frame(table(dup_ensembl$entrez_id))[,1]) #279
  tmp= unlist(by(dup_ensembl,
          dup_ensembl$entrez_id,
          function(x) x[which(!is.na(x$gene_id)),"ensembl_gene_id"]))
  Brainspan_RNAseq$entrez_id<-Brainspan_RNAseq_row$entrez_id
  tmp_table<-Brainspan_RNAseq[which(rownames(Brainspan_RNAseq) %in% tmp),1:length(Brainspan_RNAseq[1,])-1]
  tmp_table_entrez<-Brainspan_RNAseq[which(rownames(Brainspan_RNAseq) %in% tmp),length(Brainspan_RNAseq[1,])]
  tmp= as.character(as.data.frame(table(by(tmp_table,
          tmp_table_entrez,
          function(x) rownames(x)[which.max(rowMeans(x))])))$Var1)
  retain_ensembl<-c(retain_ensembl,tmp)
  dup_collect<-unlist(dup_ensembl[which(dup_ensembl$ensembl_gene_id %in% tmp),"entrez_id"])
  dup_ensembl<-dup_ensembl[!(dup_ensembl$entrez_id %in% dup_collect),]
  dup_ensembl_table<-Brainspan_RNAseq[which(Brainspan_RNAseq$entrez_id %in% dup_ensembl$entrez_id),1:ncol(Brainspan_RNAseq)-1]
  tmp= as.character(as.data.frame(table(by(dup_ensembl_table,
          dup_ensembl$entrez_id,
          function(x) rownames(x)[which.max(rowMeans(x))])))$Var1)
  retain_ensembl<-c(retain_ensembl,tmp)
  return(retain_ensembl)
}

# check variance of expression values per gene

screen_table_cluster_var<-function(table_cluster){
	variance<-c()
	for(i in 1:length(table_cluster[,1])){
		tmp<-var(as.numeric(table_cluster[i,]))
		variance<-c(variance,tmp)
	}
	table_cluster$variance<-variance
	return(table_cluster)
}

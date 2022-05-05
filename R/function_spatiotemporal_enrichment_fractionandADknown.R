# format transformation

## data format processing: permutation test

formatchange_ADknown_padjust<-function(df){
  df_c<-as.data.frame(matrix(unlist(strsplit(df[,1],split = '_')),ncol=2,byrow=T))
  colnames(df_c)<-c("cluster","win")
  df_c$p_value<-df[,2]
  df_c$p_adjust<-p.adjust(df_c$p_value,method = "fdr")
  df_c<-df_c[order(df_c$cluster),]
  if("TRUE" %in% unlist(lapply(df_c$p_adjust, function(x){x<0.1}))){
    print("exist p.adjust < 0.1")
  }
  return(df_c)
}

formatchange_fraction_padjust<-function(df){
  rownames(df)<-df[,1];df<-df[,-c(1)]
  colnames(df)<-paste0("win",c(1:11))
  df<-df[order(rownames(df)),]
  df<-as.data.frame(t(df))
  df_padj<-data.frame()
  just<-c()
  for(i in 1:5){
    if(i != 5){
      tmp<-p.adjust(df[,i],method = "fdr")
    }
    else{
      tmp<-c("NA",p.adjust(df[2:nrow(df),i],method = "fdr"))
    }
    options(warn=-1)
    tmp_df<-data.frame(rep(paste0("cluster",i),11));colnames(tmp_df)<-"cluster"
    tmp_df$win<-paste0("win",1:11)
    tmp_df$p_value<-df[,i]
    tmp_df$p_adjust<-as.numeric(tmp)
    df_padj<-rbind(df_padj,tmp_df)
    just_tmp<-unlist(lapply(tmp, function(x){tmp<0.1}))
    just<-c(just,just_tmp)
  }
  if("TRUE" %in% just){
    print("exist p.adjust < 0.1")
  }
  return(df_padj)
}

## data format processing: fraction

formatchange_fraction_calculation<-function(df){
  rownames(df)<-df[,1];df<-df[,-c(1)]
  df<-df[order(rownames(df)),]
  df<-as.data.frame(t(df))
  df_c<-data.frame()
  for(i in 1:5){
    tmp<-df[,i]
    tmp_df<-as.data.frame(tmp)
    colnames(tmp_df)<-"fraction"
    tmp_df$cluster<-rep(paste0("cluster",as.character(i)),length(tmp_df$fraction))
    tmp_df$win<-paste0("win",1:11)
    df_c<-rbind(df_c,tmp_df)
  }
  return(df_c)
}

# draw line plot

## permutation test results

draw_lineplot<-function(df,geneset.name){
  df$win <- factor(df$win,levels=paste0("win",1:11))
  y.max=max(-log10(df$p_adjust))+0.1
  ggplot(df, aes(x = win, y = -log10(p_adjust), colour = cluster)) +
    geom_point(size = 3) +
    geom_line(aes(group = cluster), size = 1) + # draw the lines according to their groups
    scale_colour_manual(values = c("#53099a","#2B66AE", "#9FCAE9", "#EA492C", "#FFC4BC")) +
    scale_y_continuous(expand = c(0.01, 0.05),limits = c(0,y.max)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 12),
      axis.ticks.length = unit(0.3, "cm"),
      axis.ticks = element_line(size = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey"),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(hjust = 0.5),
    ) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(title = geneset.name,x = "win",y = "-log10(p-adjust)",colour = "") +
    geom_hline(yintercept = -log10(0.1),linetype='dashed',colour = '#7c7c7c',size = 1)
}

## fraction calculation results

draw_lineplot_fraction_calculation<-function(df.fraction,df.p_adjust,geneset.name){
  df<-merge(df.fraction,df.p_adjust,by = c("cluster","win"))
  df$win <- factor(df$win,levels=paste0("win",1:11))
  df$sig_status<-unlist(lapply(df$p_adjust,function(x){if(!(x %in% NA)){if(x<0.1){"*"}else{""}}else{""}}))
  df$sig_color<-unlist(lapply(df$cluster,function(x){if(x=="cluster1"){"#53099a"}else if(x=="cluster2"){"#2B66AE"}else if(x=="cluster3"){"#9FCAE9"}else if(x=="cluster4"){"#EA492C"}else if(x=="cluster5"){"#FFC4BC"}}))
  y.max=max(df$fraction)+0.01
  ggplot(df, aes(x = win, y = fraction, colour = cluster)) +
    geom_point(size = 3) +
    geom_line(aes(group = cluster), size = 1) + # draw the lines according to their groups
    scale_colour_manual(values = c("#53099a","#2B66AE", "#9FCAE9", "#EA492C", "#FFC4BC")) +
    ylim(0,y.max) + # notice: here use "ylim" instead of "limits" in the "scale_y_continuous" function is more effective when in combination with "scales::percent"
    scale_y_continuous(expand = c(0.01, 0.05),labels = scales::percent) + #limits = c(0.0,y.max),
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 12),
      axis.ticks.length = unit(0.3, "cm"),
      axis.ticks = element_line(size = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey"),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(hjust = 0.5),
    ) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(title = geneset.name,x = "win",y = "Fraction of co-expressed interacting pairs",colour = "") +
    geom_text(aes(y = fraction,label=df$sig_status,vjust=0.1), size=8,face="bold",color=df$sig_color,position = position_dodge(width=0.00),check_overlap = FALSE)
}

library(openxlsx)
load("./spatiotemporal_enrichment_plots_varscreen.RData")

for(i in 1:5){
	for(j in 1:sum_windows){
		table_cluster<-get(paste("table_cluster",i,"_win",j,sep=""))
		file_name=paste("data_table_cluster/","table_cluster",i,"_win",j,".txt",sep = "")
		write.table(table_cluster,file_name,sep = "\t",col.names = T,row.names = T,quote=F)
	}
}

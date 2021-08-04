import os, sys, xlrd, xlwt
import numpy as np
import pandas as pd
import argparse
from scipy.stats import kendalltau, pearsonr, spearmanr
from multiprocessing import Pool

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--cluster", type=int, default=1)
parser.add_argument("--window", type=int, default=11)
parser.add_argument("--correlation-method", type=str, default="pearson",help='Input correlation method for calculation: pearson or spearman')
args = parser.parse_args()

os.chdir("~/spatial_temporal_analysis/")

geneset=pd.read_table("~/spatial_temporal_analysis/list_TRAPD_caseonly_mutated_an_burdengenes_entrezid.txt",header=None)
geneset=list(geneset.iloc[:, 0])

def pearsonr_pval(x, y):
    return pearsonr(x, y)[1]

def spearmanr_pval(x, y):
    return spearmanr(x, y)[1]

def window_cortest(window_j,methods):
    j=window_j
    #import expression data
    file_name="~/spatial_temporal_analysis/data_table_cluster/"+"table_cluster"+str(args.cluster)+"_win"+str(j)+".txt"
    tmp_data = pd.read_table(file_name, sep='\t')
    t_tmp_data=pd.DataFrame(tmp_data.values.T,index=tmp_data.columns,columns=tmp_data.index)
    tmp_cor=t_tmp_data.corr(method=methods)
    np_tmp_cor=tmp_cor.mask(np.triu(np.ones(tmp_cor.shape)).astype(bool)).stack()
    df_np_tmp_cor=pd.DataFrame(np_tmp_cor)
    df_np_tmp_cor_reset=df_np_tmp_cor.reset_index(level=[0,1])
    df_np_tmp_cor_reset.columns=["gene1","gene2","cor"]
    df_np_tmp_cor_reset=df_np_tmp_cor_reset[(df_np_tmp_cor_reset['gene1'].isin(geneset)) | (df_np_tmp_cor_reset['gene2'].isin(geneset))]
    tmp_spatiotemp=pd.DataFrame(["cluster"+str(args.cluster)+"_win"+str(j)]*len(df_np_tmp_cor_reset.iloc[:, 0]),columns=["spatio-temperal"])
    df_np_tmp_cor_reset.reset_index(drop=True, inplace=True)
    tmp_spatiotemp.reset_index(drop=True, inplace=True)
    df_np_tmp_cor_reset = pd.concat([df_np_tmp_cor_reset,tmp_spatiotemp],axis=1)
    return(df_np_tmp_cor_reset)

def window_cortest_pval(window_j,methods):
    j=window_j
    #import expression data
    file_name="~/spatial_temporal_analysis/data_table_cluster/"+"table_cluster"+str(args.cluster)+"_win"+str(j)+".txt"
    tmp_data = pd.read_table(file_name, sep='\t')
    t_tmp_data=pd.DataFrame(tmp_data.values.T,index=tmp_data.columns,columns=tmp_data.index)
    if methods=="pearson":
        tmp_corpval=t_tmp_data.corr(method=pearsonr_pval)
    elif methods=="spearman":
        tmp_corpval=t_tmp_data.corr(method=spearmanr_pval)
    np_tmp_corpval=tmp_corpval.mask(np.triu(np.ones(tmp_corpval.shape)).astype(bool)).stack()
    df_np_tmp_corpval=pd.DataFrame(np_tmp_corpval)
    df_np_tmp_corpval_reset=df_np_tmp_corpval.reset_index(level=[0,1])
    df_np_tmp_corpval_reset.columns=["gene1","gene2","p_value"]
    df_np_tmp_corpval_reset=df_np_tmp_corpval_reset[(df_np_tmp_corpval_reset['gene1'].isin(geneset)) | (df_np_tmp_corpval_reset['gene2'].isin(geneset))]
    tmp_spatiotemp=pd.DataFrame(["cluster"+str(args.cluster)+"_win"+str(j)]*len(df_np_tmp_corpval_reset.iloc[:, 0]),columns=["spatio-temperal"])
    df_np_tmp_corpval_reset.reset_index(drop=True, inplace=True)
    tmp_spatiotemp.reset_index(drop=True, inplace=True)
    df_np_tmp_corpval_reset = pd.concat([df_np_tmp_corpval_reset,tmp_spatiotemp],axis=1)
    return(df_np_tmp_corpval_reset)

if __name__ == "__main__":
    # parallel running the correlation test
    ## Correlation coefficient
    pool = Pool()
    jobs = []
    for window_j in range(1,args.window+1):
        jobs.append(pool.apply_async(window_cortest, args=(window_j,args.correlation_method,)))
    pool.close()
    pool.join()
    cor_data_test = pd.DataFrame()
    for i in jobs:
        tmp_dataframe = i.get()
        cor_data_test = pd.concat([cor_data_test, tmp_dataframe], axis=0)
    # return the coefficient matrix
    filename_cor_data_test="~/spatial_temporal_analysis/padj_correlation_adjust_caseonly/"+"caseonlyburden_coexp_cor_"+"cluster"+str(args.cluster)+"_"+args.correlation_method+".txt"
    cor_data_test.to_csv(filename_cor_data_test,sep='\t',index=False)
    ## Correlation p-values
    #del cor_data_test,tmp_dataframe
    pool = Pool()
    jobs = []
    for window_j in range(1,args.window+1):
        jobs.append(pool.apply_async(window_cortest_pval, args=(window_j,args.correlation_method,)))
    pool.close()
    pool.join()
    cor_data_test_pval = pd.DataFrame()
    for i in jobs:
        tmp_dataframe_pval = i.get()
        cor_data_test_pval = pd.concat([cor_data_test_pval, tmp_dataframe_pval], axis=0)
    ## merge cor and pval parameters for the same gene pairs
    # return the intact matrix
    cor_data_test_all=pd.merge(cor_data_test,cor_data_test_pval,on=["gene1","gene2","spatio-temperal"],how='outer')
    filename_cor_data_test_all="~/spatial_temporal_analysis/padj_correlation_adjust_caseonly/"+"caseonlyburden_coexp_all_"+"cluster"+str(args.cluster)+"_"+args.correlation_method+".txt"
    cor_data_test_all.to_csv(filename_cor_data_test_all,sep='\t',index=False)

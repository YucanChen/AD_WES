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
parser.add_argument("--geneset", type=str, help='Input the file path to the geneset.txt for correlation network building')
args = parser.parse_args()

os.chdir("~/spatial_temporal_analysis/")

geneset=pd.read_table(args.geneset, header=None)
geneset=list(geneset.iloc[:, 0])

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

if __name__ == "__main__":
    # parallel running the correlation test
    ## Correlation coefficient
    pool = Pool()
    jobs = []
    for window_j in range(1,args.window+1):# Note: cluster5 range(2,args.window+1)
        jobs.append(pool.apply_async(window_cortest, args=(window_j,args.correlation_method,)))
    pool.close()
    pool.join()
    cor_data_test = pd.DataFrame()
    for i in jobs:
        tmp_dataframe = i.get()
        cor_data_test = pd.concat([cor_data_test, tmp_dataframe], axis=0)
    # return the coefficient matrix
    filename_cor_data_test="~/spatial_temporal_analysis/padj_correlation_adjust/"+"known_ADgenes_v2_coexp_subcorr_"+"cluster"+str(args.cluster)+"_"+args.correlation_method+".txt"
    cor_data_test.to_csv(filename_cor_data_test,sep='\t',index=False)

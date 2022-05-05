import os, sys, xlrd, xlwt
import pandas as pd
import numpy as np
import argparse

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--cluster', type=int, help='Input the number of the cluster')
parser.add_argument("--window", type=int, help='Input the number of the window')
parser.add_argument('--geneset', type=str, help='Input geneset for permutation test in AD known geneset coexpression network')
parser.add_argument('--correlation-method', type=str, help='the correlation method used for network construction: pearson or spearman')
parser.add_argument('--correlation-threshold', type=float, help='Input the correlation threshold for PPI pairs')
args = parser.parse_args()

os.chdir('./results/')

if __name__ == '__main__':
  geneset_list=pd.read_table(args.geneset,header=None)
  geneset_list=list(geneset_list.iloc[:, 0])
  file_name='./results/net_thresholdcor_PPI/'+'net_coexp_threshold'+str(args.correlation_threshold)+'_'+args.correlation_method+'_cluster'+str(args.cluster)+'.txt'
  tmp_data = pd.read_table(file_name, sep='\t')
  all_geneset_data = tmp_data[(tmp_data['gene1'].isin(geneset_list)) | (tmp_data['gene2'].isin(geneset_list))]
  spatio_temperal = "cluster"+str(args.cluster)+"_win"+str(args.window)
  tmp_data_temp = all_geneset_data[all_geneset_data['spatio-temperal'] == spatio_temperal]
  #output the screened result
  filename_screen_data='./results/net_thresholdcor_PPI_specificsets/'+'geneset_'+'net_coexp_threshold'+str(args.correlation_threshold)+'_'+args.correlation_method+'_cluster'+str(args.cluster)+'_win'+str(args.window)+'.txt'
  tmp_data_temp.to_csv(filename_screen_data,sep='\t',index=False)

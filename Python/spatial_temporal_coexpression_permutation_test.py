import os, sys, xlrd, xlwt
import pandas as pd
import numpy as np
import argparse
import random
from multiprocessing import Pool

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--cluster", type=int, default=1)
parser.add_argument("--window", type=int, default=1)
parser.add_argument("--permutation-time", type=int, default=1000)
args = parser.parse_args()

file_name="~/spatial_temporal_analysis/padj_correlation_adjust/"+"known_ADgenes_coexp_subcorr_cluster"+str(args.cluster)+".txt"

#import coexpression data
data = pd.read_table(file_name, sep='\t')
known_ADlist = pd.read_table("AD_knownsig_geneset_v1.txt", header=None)
known_ADlist = list(known_ADlist.iloc[:, 0])
burdengenes_AD_list=pd.read_table("~/spatial_temporal_analysis/burdengenes_AD_list.txt",header=None)
burdengenes_AD_list=list(burdengenes_AD_list.iloc[:, 0])
# calculate number of expressed burden genes in each network
spatiotemperal = "cluster" + str(args.cluster) + "_" + "win" + str(args.window)
data_test = data[data['spatio-temperal'] == spatiotemperal]
#calculate number of burdengenes_AD
burdengenes_AD_num = 0
cor_data_test = pd.DataFrame
tmp_dataframe = pd.DataFrame
for i in range(len(burdengenes_AD_list)):
    if (burdengenes_AD_list[i] in list(data_test['gene1'])) or (burdengenes_AD_list[i] in list(data_test['gene2'])):
        burdengenes_AD_num += 1
        tmp_dataframe1 = data_test[(data_test['gene1'] == burdengenes_AD_list[i])]
        tmp_dataframe2 = data_test[(data_test['gene2'] == burdengenes_AD_list[i])]
        cor_data_test = pd.concat([tmp_dataframe1, tmp_dataframe2], axis=0)
# output median of burden genes
burdengenes_AD_median = np.median(abs(cor_data_test["cor"]))
# generate random seed genes which is the same number as the expressed genes in each network
random_gene_list = list(set(list(data_test['gene1']) + list(data_test['gene2'])))
# main_list = [item for item in random_gene_list if item not in list(set(known_ADlist+burdengenes_AD_list))]
random_gene_list = list(set(random_gene_list) - set(known_ADlist + burdengenes_AD_list))

def permutation_test(permutation_time):
    #for j in range(0,args.permutation_time):
    j=permutation_time
    random_seeds = random.sample(random_gene_list, burdengenes_AD_num)
    cor_data_test=pd.DataFrame
    tmp_dataframe=pd.DataFrame
    for i in range(0,burdengenes_AD_num):
        if (random_seeds[i] in list(data_test['gene1'])) or (random_seeds[i] in list(data_test['gene2'])):
            tmp_dataframe1=data_test[(data_test['gene1']==random_seeds[i])]
            tmp_dataframe2 = data_test[(data_test['gene2'] == random_seeds[i])]
            cor_data_test=pd.concat([tmp_dataframe1,tmp_dataframe2],axis=0)
    tmp_median=np.median(abs(cor_data_test["cor"]))
    return(tmp_median)

if __name__ == "__main__":
  # calculate median value
  # iterate multiple times and obtain the distribution of medians
  pool = Pool()
  jobs = []
  for permutation_time in range(args.permutation_time):
    jobs.append(pool.apply_async(permutation_test, args=(permutation_time,)))
  pool.close()
  pool.join()
  random_median_list = []
  for i in jobs:
    random_median_list.append(i.get())

#significance test
p_value=(len([item for item in random_median_list if item>burdengenes_AD_median])+1)/(args.permutation_time+1)
#return the significant results
print("cluster"+str(args.cluster)+"_"+"win"+str(args.window)+'\t'+str(p_value))

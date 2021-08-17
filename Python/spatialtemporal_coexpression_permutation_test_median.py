import os, sys, xlrd, xlwt
import pandas as pd
import numpy as np
import argparse
from multiprocessing import Pool

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--cluster", type=int, default=1)
parser.add_argument("--window", type=int, default=1)
parser.add_argument("--permutation-time", type=int, default=1000)
parser.add_argument("--geneset", type=str, help='Input geneset for permutation test in AD known geneset coexpression network')
parser.add_argument("--fitlength-distributiongenes", type=str, help='Input genes fitted the length distribution as the specific geneset')
parser.add_argument("--correlation-method", type=str, help='the correlation method used for network construction: pearson or spearman')
args = parser.parse_args()

os.chdir("~/spatial_temporal_analysis/")
file_name="~/spatial_temporal_analysis/padj_correlation_adjust/"+"known_ADgenes_v2_coexp_subcorr_cluster"+str(args.cluster)+"_"+args.correlation_method+".txt"

#import coexpression data
data = pd.read_table(file_name, sep='\t')
known_ADlist = pd.read_table("AD_knownsig_geneset_v2.txt", header=None)
known_ADlist = list(known_ADlist.iloc[:, 0])
geneset_list=pd.read_table(args.geneset,header=None)
geneset_list=list(geneset_list.iloc[:, 0])
# calculate number of expressed genes from geneset in each network
spatiotemperal = "cluster" + str(args.cluster) + "_" + "win" + str(args.window)
data_test = data[data['spatio-temperal'] == spatiotemperal]

# calculate number of genes from the geneset
geneset_num = 0
cor_data_test = pd.DataFrame()
tmp_dataframe = pd.DataFrame()
for i in range(len(geneset_list)):
    if (geneset_list[i] in list(data_test['gene1'])) or (geneset_list[i] in list(data_test['gene2'])):
        geneset_num += 1
        tmp_dataframe1 = data_test[(data_test['gene1'] == geneset_list[i])]
        tmp_dataframe2 = data_test[(data_test['gene2'] == geneset_list[i])]
        tmp_dataframe=pd.concat([tmp_dataframe1,tmp_dataframe2],axis=0)
        cor_data_test=pd.concat([cor_data_test,tmp_dataframe],axis=0)
# output median of genes from the geneset
geneset_median = np.median(abs(cor_data_test["cor"]))

# generate random seed genes which is the same number as the expressed genes from the geneset in each network
fitlength_distribution_genes=pd.read_table(args.fitlength_distributiongenes)
random_gene_list = list(set(list(data_test['gene1']) + list(data_test['gene2'])))
random_gene_list = list(set(random_gene_list).intersection(set(fitlength_distribution_genes['entrezgene_id'])))
random_gene_list = list(set(random_gene_list) - set(known_ADlist + geneset_list))

def permutation_test(permutation_time):
    j=permutation_time
    random_seeds = list(np.random.choice(random_gene_list, geneset_num, replace=False))
    cor_data_test=pd.DataFrame()
    tmp_dataframe=pd.DataFrame()
    for i in range(0,geneset_num):
        if (random_seeds[i] in list(data_test['gene1'])) or (random_seeds[i] in list(data_test['gene2'])):
            tmp_dataframe1 = data_test[(data_test['gene1']==random_seeds[i])]
            tmp_dataframe2 = data_test[(data_test['gene2'] == random_seeds[i])]
            tmp_dataframe=pd.concat([tmp_dataframe1,tmp_dataframe2],axis=0)
            cor_data_test=pd.concat([cor_data_test,tmp_dataframe],axis=0)
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
p_value=(len([item for item in random_median_list if item > geneset_median])+1)/(args.permutation_time+1)
#return the significant results
print("cluster"+str(args.cluster)+"_"+"win"+str(args.window)+'\t'+str(p_value))

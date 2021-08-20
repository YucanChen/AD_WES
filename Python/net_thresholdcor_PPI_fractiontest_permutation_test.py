import os, sys, xlrd, xlwt
import pandas as pd
import numpy as np
import argparse
from multiprocessing import Pool

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--permutation-time', type=int, default=1000)
parser.add_argument('--cluster', type=int, help='Input the number of the cluster')
parser.add_argument('--geneset', type=str, help='Input geneset for permutation test in AD known geneset coexpression network')
parser.add_argument('--fitlength-distributiongenes', type=str, help='Input genes fitted the length distribution as the specific geneset')
parser.add_argument('--correlation-method', type=str, help='the correlation method used for network construction: pearson or spearman')
parser.add_argument('--correlation-threshold', type=float, help='Input the correlation threshold for PPI pairs')
args = parser.parse_args()

os.chdir('~/spatial_temporal_analysis/')
geneset_list=pd.read_table(args.geneset,header=None)
geneset_list=list(geneset_list.iloc[:, 0])
geneset_num=len(geneset_list)

# calculate fraction of co-expression pairs of geneset according to different windows in the same cluster
file_name='~/spatial_temporal_analysis/net_thresholdcor_PPI/'+'net_coexp_threshold'+str(args.correlation_threshold)+'_'+args.correlation_method+'_cluster'+str(args.cluster)+'.txt'
tmp_data = pd.read_table(file_name, sep='\t')
all_geneset_data = tmp_data[(tmp_data['gene1'].isin(geneset_list)) | (tmp_data['gene2'].isin(geneset_list))]
num_all_cluster_pairs = len(all_geneset_data['cor'])
if args.cluster != 5:
  start = 1
  geneset_fraction_win_list = []
else:
  start = 2
  geneset_fraction_win_list = ['NA']
for j in range(start,12):
  spatio_temperal = 'cluster'+str(args.cluster)+'_win'+str(j)
  tmp_data_temp = all_geneset_data[all_geneset_data['spatio-temperal'] == spatio_temperal]
  num_tmp_data_temp = len(tmp_data_temp['cor'])
  geneset_fraction_win_list.append(num_tmp_data_temp/num_all_cluster_pairs)

# generate random seed genes
fitlength_distribution_genes=pd.read_table(args.fitlength_distributiongenes)
random_gene_list = list(set(fitlength_distribution_genes['entrezgene_id']))
random_gene_list = list(set(random_gene_list) - set(geneset_list))

def permutation_test(permutation_time):
  m=permutation_time
  random_seeds = list(np.random.choice(random_gene_list, geneset_num, replace=False))
  all_random_data = tmp_data[(tmp_data['gene1'].isin(random_seeds)) | (tmp_data['gene2'].isin(random_seeds))]
  num_all_cluster_pairs = len(all_random_data['cor'])
  if args.cluster == 5:
    random_fraction_win_list = ['NA']
  else:
    random_fraction_win_list = []
  for win in range(start,12):
    spatio_temperal = 'cluster'+str(args.cluster)+'_win'+str(win)
    tmp_data_temp = all_random_data[all_random_data['spatio-temperal'] == spatio_temperal]
    num_tmp_data_temp = len(tmp_data_temp['cor'])
    random_fraction_win_list.append(num_tmp_data_temp/num_all_cluster_pairs)
  return(random_fraction_win_list)

if __name__ == '__main__':
  cpu_num=os.cpu_count()
  pool = Pool(cpu_num)
  jobs = []
  jobs_res = []
  for permutation_time in range(args.permutation_time):
    jobs.append(pool.apply_async(permutation_test, args=(permutation_time,)))
  pool.close()
  pool.join()
  for permu in jobs:
    jobs_res.append(permu.get())
  jobs_res = pd.DataFrame(jobs_res)
  jobs_res.columns = ['win'+ str(win) for win in range(1,12)]
  if args.cluster == 5:
    cluster_res = ['cluster'+str(args.cluster),'NA']
  else:
    cluster_res = ['cluster'+str(args.cluster)]
  for j in range(start,12):
    permu_win_list = list(jobs_res['win'+str(j)])
    #significance test
    permu_win_pvalue=(len([item for item in permu_win_list if item > geneset_fraction_win_list[j-1]])+1)/(args.permutation_time+1)
    cluster_res.append(permu_win_pvalue)
  
  #output the permutation results
  outfile='~/spatial_temporal_analysis/net_thresholdcor_PPI_fractiontest/fraction_permutationtest_'+'net_coexp_threshold'+str(args.correlation_threshold)+'_'+args.correlation_method+'_'+str(args.permutation_time)+'.txt'
  res_out=str(cluster_res).replace('[','').replace(']','').replace("'",'').replace(', ','\t')
  f=open(outfile,'a')
  f.write(res_out)
  f.write('\n')
  f.close()

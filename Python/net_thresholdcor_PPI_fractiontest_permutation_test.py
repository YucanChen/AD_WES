import os, sys, xlrd, xlwt
import pandas as pd
import numpy as np
import argparse
from multiprocessing import Pool

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--permutation-time", type=int, default=1000)
parser.add_argument("--geneset", type=str, help='Input geneset for permutation test in AD known geneset coexpression network')
parser.add_argument("--fitlength-distributiongenes", type=str, help='Input genes fitted the length distribution as the specific geneset')
parser.add_argument("--correlation-method", type=str, help='the correlation method used for network construction: pearson or spearman')
parser.add_argument("--correlation-threshold", type=float, help='Input the correlation threshold for PPI pairs')
args = parser.parse_args()

os.chdir("~/spatial_temporal_analysis/")
geneset_list=pd.read_table(args.geneset,header=None)
geneset_list=list(geneset_list.iloc[:, 0])
geneset_num=len(geneset_list)

# calculate fraction of co-expression pairs of geneset according to different windows in the same cluster
geneset_fraction_cluster_list = []
for i in range(1,6):
	file_name="~/spatial_temporal_analysis/net_thresholdcor_PPI/"+"net_coexp_threshold"+str(args.correlation_threshold)+"_"+args.correlation_method+"_cluster"+str(i)+".txt"
	tmp_data = pd.read_table(file_name, sep='\t')
	all_geneset_data = tmp_data[(tmp_data['gene1'].isin(geneset_list)) | (tmp_data['gene2'].isin(geneset_list))]
	num_all_cluster_pairs = len(all_geneset_data['cor'])
	if i != 5:
		s=1
		geneset_fraction_tmp_list = []
	else:
		s=2
		geneset_fraction_tmp_list = ["NA"]
	tmp_data_temp = pd.DataFrame()
	for j in range(s,12):
		spatio_temperal = "cluster"+str(i)+"_win"+str(j)
		tmp_data_temp = all_geneset_data[all_geneset_data['spatio-temperal'] == spatio_temperal]
		num_tmp_data_temp = len(tmp_data_temp['cor'])
		geneset_fraction_tmp_list.append(num_tmp_data_temp/num_all_cluster_pairs)
	geneset_fraction_cluster_list.append(geneset_fraction_tmp_list)
geneset_fraction_cluster_df = pd.DataFrame(geneset_fraction_cluster_list)
geneset_fraction_cluster_df.columns = ["win1","win2","win3","win4","win5","win6","win7","win8","win9","win10","win11"]
geneset_fraction_cluster_df.index = ["cluster1","cluster2","cluster3","cluster4","cluster5"]

# generate random seed genes
fitlength_distribution_genes=pd.read_table(args.fitlength_distributiongenes)
random_gene_list = list(set(fitlength_distribution_genes['entrezgene_id']))
random_gene_list = list(set(random_gene_list) - set(geneset_list))

def permutation_test(permutation_time):
  j=permutation_time
  random_seeds = list(np.random.choice(random_gene_list, geneset_num, replace=False))
  random_fraction_cluster_list = []
  for i in range(1,6):
    file_name="~/spatial_temporal_analysis/net_thresholdcor_PPI/net_coexp_threshold"+str(args.correlation_threshold)+"_"+args.correlation_method+"_cluster"+str(i)+".txt"
    tmp_data = pd.read_table(file_name, sep='\t')
    all_random_data = tmp_data[(tmp_data['gene1'].isin(random_seeds)) | (tmp_data['gene2'].isin(random_seeds))]
    num_all_cluster_pairs = len(all_random_data['cor'])
    if i != 5:
      s=1
      random_fraction_tmp_list = []
    else:
      s=2
      random_fraction_tmp_list = ["NA"]
    tmp_data_temp = pd.DataFrame()
    for j in range(s,12):
      spatio_temperal = "cluster"+str(i)+"_win"+str(j)
      tmp_data_temp = all_random_data[all_random_data['spatio-temperal'] == spatio_temperal]
      num_tmp_data_temp = len(tmp_data_temp['cor'])
      random_fraction_tmp_list.append(num_tmp_data_temp/num_all_cluster_pairs)
    random_fraction_cluster_list.append(random_fraction_tmp_list)
  random_fraction_cluster_df = pd.DataFrame(random_fraction_cluster_list)
  random_fraction_cluster_df.columns = ["win1","win2","win3","win4","win5","win6","win7","win8","win9","win10","win11"]
  random_fraction_cluster_df.index = ["cluster1","cluster2","cluster3","cluster4","cluster5"]
  return(random_fraction_cluster_df)

if __name__ == "__main__":
	cpu_num=os.cpu_count()
	pool = Pool(cpu_num)
	jobs = []
	for permutation_time in range(args.permutation_time):
		jobs.append(pool.apply_async(permutation_test, args=(permutation_time,)))
	pool.close()
	pool.join()
	permutation_res = []
	for k in range(1,6):
		if k != 5:
			s = 1
			cluster_res = []
		else:
			s = 2
			cluster_res = ["NA"]
		for h in range(s,12):
			permu_tmp_list = []
			for i in jobs:
				tmp_df = i.get()
				permu_tmp_list.append(tmp_df["win"+str(h)]["cluster"+str(k)])
				#significance test
			permu_tmp_pvalue=(len([item for item in permu_tmp_list if item > geneset_fraction_cluster_df["win"+str(h)]["cluster"+str(k)]])+1)/(args.permutation_time+1)
			cluster_res.append(permu_tmp_pvalue)
		permutation_res.append(cluster_res)
	permutation_res_df = pd.DataFrame(permutation_res)
	permutation_res_df.columns = ["win1","win2","win3","win4","win5","win6","win7","win8","win9","win10","win11"]
	permutation_res_df.index = ["cluster1","cluster2","cluster3","cluster4","cluster5"]

	#output the permutation results
	outfile="~/spatial_temporal_analysis/net_thresholdcor_PPI_fractiontest/fraction_permutationtest_"+"net_coexp_threshold"+str(args.correlation_threshold)+"_"+args.correlation_method+"_"+str(args.permutation_time)+".csv"
	permutation_res_df.to_csv(outfile)

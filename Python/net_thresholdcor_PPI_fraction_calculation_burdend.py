import os, sys, xlrd, xlwt
import pandas as pd
import argparse

#set up parameters
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--geneset", type=str, help='Input geneset for permutation test in AD known geneset coexpression network')
parser.add_argument("--correlation-method", type=str, help='the correlation method used for network construction: pearson or spearman')
parser.add_argument("--correlation-threshold", type=float, help='Input the correlation threshold for PPI pairs')
args = parser.parse_args()

geneset_list=pd.read_table(args.geneset,header=None)
geneset_list=list(geneset_list.iloc[:, 0])

if __name__ == "__main__":
	# calculate fraction of co-expression pairs of geneset according to different windows in the same cluster
	geneset_fraction_cluster_list = []
	for i in range(1,6):
		file_name="./results/net_thresholdcor_PPI/"+"net_coexp_threshold"+str(args.correlation_threshold)+"_"+args.correlation_method+"_cluster"+str(i)+".txt"
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
	geneset_fraction_cluster_df.columns = ['win'+ str(win) for win in range(1,12)]
	geneset_fraction_cluster_df.index = ['cluster' + str(cluster) for cluster in range(1,6)]

	#output the permutation results
	outfile="./results/net_thresholdcor_PPI_fractiontest/fraction_calculation_burdend_"+"net_coexp_threshold"+str(args.correlation_threshold)+"_"+args.correlation_method+".csv"
	geneset_fraction_cluster_df.to_csv(outfile)

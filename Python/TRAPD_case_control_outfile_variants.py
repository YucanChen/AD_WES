# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:59 2021
@author: CYC
AD_exome - Assignment.5-0 conditional analysis-split the outfile variants for statistical analysis
"""
import os, sys, xlrd, xlwt
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("--i", help="input the name of the file that you want to change format")
parser.add_argument("--o", help="output the file being processed")
args = parser.parse_args()

# input TRAPD outfile data:
file_outfile=pd.read_table(args.i)

file_outfile_snps=pd.DataFrame()
for i in range(len(file_outfile["SNPS"])):
    list_snp=file_outfile.loc[i,"SNPS"].split(',')
    list_gene=list(np.repeat(file_outfile.loc[i,"#GENE"],len(list_snp)))
    table_i=pd.DataFrame({'SNPS':list_snp,'GENE':list_gene})
    file_outfile_snps=pd.concat([file_outfile_snps,table_i],axis=0)

file_outfile_snps.to_csv(args.o,index=False)

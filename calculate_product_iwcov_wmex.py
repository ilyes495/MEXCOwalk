import os.path
import pandas as pd
import csv

filepath_iwcov = "../hint/out/evaluation_tab/iwavg_cov.txt"
filepath_wmex = "../hint/out/evaluation_tab/wavg_mutex.txt"
filepath_iwcov_wmex = "../hint/out/evaluation_tab/iwcov_wmex.txt"

top_genes_to_read = 1000


iwcov_df = pd.read_table(filepath_iwcov, nrows=top_genes_to_read/100, index_col=0)
wmex_df = pd.read_table(filepath_wmex,nrows=top_genes_to_read/100, index_col=0)
mult = iwcov_df*wmex_df
print iwcov_df
print wmex_df
print mult

mult.to_csv(filepath_iwcov_wmex, sep='\t')
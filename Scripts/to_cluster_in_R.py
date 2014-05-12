import os, time
import numpy as np
import cluster_methods as clust



cancers = ["brca","coad","gbm","kirc","luad","lusc","ov","ucec"]
source_path = os.getcwd()+"/Output/"
source_path = "/Users/rathikannan/Documents/hm_27k_profile/Output/"

for cancer in cancers:

	# Read probe dict from txt
	probe_dict = clust.read_probe_dict(source_path+cancer+"/probe_dict_dis.txt")	
	
	# Convert dictionary to matrix
	# Rows are samples
	# Columns are probe IDs
	matrix_1 = clust.dict_to_matrix_1(dict)
	
	# Calculate hamming distance between samples
	hamming_1 = clust.hamming_full_csv(matrix_1,out_path+cancer+"/probe_dis_hamming_samples.txt")

 

		

	
	
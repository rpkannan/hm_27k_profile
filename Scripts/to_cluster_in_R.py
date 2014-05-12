import os, time
import numpy as np
import cluster_methods as clust



cancers = ["brca","coad","gbm","kirc","luad","lusc","ov","ucec"]
source_path = os.getcwd()+"/Output/"

for cancer in cancers:	
	dict = {}
	f_open = open(out_path+"brca/probe_dict_dis.txt","r").readlines()
	for line in f_open:
		line = line.strip().split('\t')
		if line[0] == "Samples":
			continue
		else:
			dict[line[0]] = [float(i) for i in line[1:]]
	
	# Convert dictionary to matrix
	# Rows are samples
	# Columns are probe IDs
	matrix_1 = clust.dict_to_matrix_1(dict)
	
	# Calculate hamming distance between samples
	hamming_1 = clust.hamming_full_csv(matrix_1,out_path+cancer+"/probe_dis_hamming_samples.txt")

 

		

	
	
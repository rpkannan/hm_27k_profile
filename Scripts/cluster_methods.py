import os
import cPickle as pickle
import numpy as np

############################## DATA STRUCTURE CONVERSION #################################
'''
Converts a dictionary to a matrix
Dictionary inputs: probe_dict, gene_dict
Dictionary format:
Key: barcode
Value: list of gene or probe beta values; should already be sorted in proper order
Return a numpy data matrix
Rows: indexed by barcode
Columns: indexed by probe/gene ID
'''
def dict_to_matrix_1(dict):
	barcodes = sorted(dict.keys())
	data_matrix = []
	for b in barcodes:
		data_matrix.append(dict[b])
	#convert native data array into a numpy array
	data_matrix = np.array(data_matrix)
	return data_matrix

'''
Converts a dictionary to a matrix
Dictionary inputs: probe_dict, gene_dict
Dictionary format:
Key: barcode
Value: list of gene or probe beta values; should already be sorted in proper order
Return a numpy data matrix
Rows: indexed by probe/gene ID
Columns: indexed by barcode
'''
def dict_to_matrix_2(dict):
	data_matrix = dict_to_matrix_1(dict)
	data_matrix = data_matrix.transpose()
	return data_matrix

############################## READ DICT FROM TXT ##################################
def read_probe_dict(file_name):
	dict = {}
	f_open = open(file_name,"r").readlines()
	for line in f_open:
		line = line.strip().split('\t')
		if line[0] == "Samples":
			continue
		else:
			dict[line[0]] = [float(i) for i in line[1:]]
	return dict

################################ WRITE PROBE DICT TO TXT #################################
probes_27_short = open("/Users/rathikannan/Documents/hm_27k_profile/Annotation/probes_27_short.txt","r").readline().strip().split('\t')
def probe_dict_to_txt(probe_dict,file_name):
	f_out = open(file_name,"w")
	f_line = ["Samples"] + probes_27_short
	f_out.write('\t'.join(f_line)+'\n')
	for barcode in probe_dict:
		f_out.write('\t'.join([str(i) for i in probe_dict[barcode]])+'\n')
	f_out.close()
	return

################################# FREQUENCY MEASURES #####################################

def max_freq(np_ls):
	sums = [0]*4
	sums[0] = sum(np_ls == 0)
	sums[1] = sum(np_ls == 1)
	sums[2] = sum(np_ls == 2)
	sums[3] = sum(np.isnan(np_ls))
	max_count = np.amax(sums)
	return float(max_count)/len(np_ls)

def freq(np_ls):
	sums = [0]*4
	sums[0] = sum(np_ls == 0)
	sums[1] = sum(np_ls == 1)
	sums[2] = sum(np_ls == 2)
	sums[3] = sum(np.isnan(np_ls))
	return sums

'''
Calculates a resolved methylation score.

Input: a list of discretized beta values
Output: single float score of value: 0.0,1.0,2.0,np.nan

Scoring Algorithm: 
If input list has only one value, that values is returned.
If input list has 2+ values, the score with the highest frequency of occurrence is returned. 
If there is a tie, the following hierarchy is used as the tie break: 2 > 0 > 1 > NA
'''
def resolve(ls):
	ls = np.array(ls)
	sums = [0]*4
	sums[0] = sum(ls == 0)
	sums[1] = sum(ls == 1)
	sums[2] = sum(ls == 2)
	sums[3] = sum(np.isnan(ls))
	
	max_val = np.argwhere(sums == np.amax(sums))
	
	if len(max_val) > 1:
		if 2 in max_val:
			return 2.0
		elif 0 in max_val:
			return 0.0
		elif  1 in max_val:
			return 1.0
		else:
			return np.nan
	else:
		if max_val[0,0] == 3:
			return np.nan
		else:
			return float(max_val[0,0])

################################## HAMMING DISTANCE ######################################
'''
Calculates hamming distance between two numpy 1D arrays of equal length
'''
def hamming_distance(s1, s2):
	if len(s1) != len(s2):
		raise ValueError("Undefined for sequences of unequal length")
	return (s1 != s2).mean() * len(s1)

'''
Calculates the pairwise hamming distances between the rows of the input data_matrix
Returns a full distance matrix (n x n for n = len(data_matrix))
Output is pickled to file_name
'''
def hamming_full_pickle(data_matrix,file_name):
	nrow = len(data_matrix)
	mat_dist = np.zeros(shape=(len(data_matrix),len(data_matrix)))
	for r in range(1,nrow):
		for c in range(0,r):
			mat_dist[r][c] = hamming_distance(data_matrix[r],data_matrix[c])
			mat_dist[c][r] = mat_dist[r][c]
	if file_name == "NA":
		return mat_dist
	else:
		pickle.dump(mat_dist,open(file_name,"wb"))
		return mat_dist

'''
Calculates the pairwise hamming distances between the rows of the input data_matrix
Returns a lower distance matrix (jagged matrix)
Output is pickled to file_name
'''
def hamming_lower_pickle(data_matrix,file_name):
	nrow = len(data_matrix)
	lower_dist = []
	for r in range(1,nrow):
		values = []
		for c in range(0,r):
			values.append(hamming_distance(data_matrix[r],data_matrix[c]))
		lower_dist.append(values)
	if file_name == "NA":
		return lower_dist
	else:
		pickle.dump(lower_dist,open(file_name,"wb"))
		return lower_dist

def hamming_lower_ls(data_matrix):
	nrow = len(data_matrix)
	lower_dist = []
	for r in range(1,nrow):
		for c in range(0,r):
			lower_dist.append(hamming_distance(data_matrix[r],data_matrix[c]))
	return lower_dist

'''
Calculates the pairwise hamming distances between the rows of the input data_matrix
Returns a full distance matrix (n x n for n = len(data_matrix))
Output is written to a tab delimited file
'''
def hamming_full_csv(data_matrix,file_name):
	nrow = len(data_matrix)
	mat_dist = np.zeros(shape=(len(data_matrix),len(data_matrix)))
	for r in range(1,nrow):
		for c in range(0,r):
			mat_dist[r][c] = hamming_distance(data_matrix[r],data_matrix[c])
			mat_dist[c][r] = mat_dist[r][c]
	if file_name == "NA":
		return mat_dist
	else:
		file_out = open(file_name,"w")
		for i in mat_dist:
			file_out.write("\t".join(str(x) for x in i))
			file_out.write('\n')
		return mat_dist

'''
Calculates the pairwise hamming distances between the rows of the input data_matrix
Returns a lower distance matrix (n x n for n = len(data_matrix))
Output is written to a tab delimited file
'''
def hamming_lower_csv(data_matrix, file_name):
	nrow = len(data_matrix)
	f_out = open(file_name,"w")
	for r in range(1,nrow):
		for c in range(0,r):
			f_out.write(str(hamming_distance(data_matrix[r],data_matrix[c]))+'\t')
		f_out.write('\n')
	f_out.close()
	return


	
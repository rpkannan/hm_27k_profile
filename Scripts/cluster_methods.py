import os
import cPickle as pickle
import numpy as np
import scipy.cluster.hierarchy as hier
import matplotlib.pylab as plt

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
	data_matrix = read_to_matrix_1(dict)
	data_matrix = data_matrix.transpose()
	return data_matrix

############################## READ DICT FROM TXT ##################################

def read_probe_dict(file_name)

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

############################## CLUSTER & VISUALIZE #######################################	
# Uses a reduced distance matrix to perform hierarchical clustering
# Returns a linkage matrix
def cluster(red_dist_matrix,i_method):
	linkage_matrix = hier.linkage(red_dist_matrix, method=i_method)
	return linkage_matrix

# Creates a distance labeled dendrogram from an input linkage matrix
def augmented_dendrogram(*args, **kwargs):
    ddata = hier.dendrogram(*args, **kwargs)
    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            plt.plot(x, y, 'ro')
            plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),textcoords='offset points',va='top', ha='center')
    return ddata

def save_image(file_name):
	plt.savefig(file_name)
	return


	  
	
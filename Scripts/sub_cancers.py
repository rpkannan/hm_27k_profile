import os
import numpy as np
import cluster_methods as clust


# Construct sub_cancers dictionaries
cancers = ["brca","gbm","kirc","luad"]


for cancer in cancers:
	
	print "Started: "+str(cancer)
	
	f_path = os.getcwd()+"Output/"+str(cancer)+"/"
	f_path = "/Users/rathikannan/Documents/hm_27k_profile/Output/"+str(cancer)+"/"

	cut1 = open(f_path+"cut1.csv","r").readlines()
	cut1 = [i.strip() for i in cut1]
	cut2 = open(f_path+"cut2.csv","r").readlines()
	cut2 = [i.strip() for i in cut2]

	probe_dict = clust.read_probe_dict(f_path+"probe_dict_dis.txt")
	dict1 = {}
	dict2 = {}
	for key in probe_dict:
		if key in cut1:
			dict1[key] = probe_dict[key]
		else:
			dict2[key] = probe_dict[key]
	
	clust.probe_dict_to_txt(dict1,f_path+str(cancer)+"_1.txt")
	clust.probe_dict_to_txt(dict2,f_path+str(cancer)+"_2.txt")



cancer = "ov"

print "Started: "+str(cancer)

f_path = os.getcwd()+"Output/"+str(cancer)+"/"
f_path = "/Users/rathikannan/Documents/hm_27k_profile/Output/"+str(cancer)+"/"

cut1 = open(f_path+"cut1.csv","r").readlines()
cut1 = [i.strip() for i in cut1]
cut2 = open(f_path+"cut2.csv","r").readlines()
cut2 = [i.strip() for i in cut2]
cut3 = open(f_path+"cut3.csv","r").readlines()
cut3 = [i.strip() for i in cut3]
cut4 = open(f_path+"cut4.csv","r").readlines()
cut4 = [i.strip() for i in cut4]

probe_dict = clust.read_probe_dict(f_path+"probe_dict_dis.txt")
dict1 = {}
dict2 = {}
dict3 = {}
dict4 = {}
for key in probe_dict:
	if key in cut1:
		dict1[key] = probe_dict[key]
	elif key in cut2:
		dict2[key] = probe_dict[key]
	elif key in cut3:
		dict3[key] = probe_dict[key]
	else:
		dict4[key] = probe_dict[key]

clust.probe_dict_to_txt(dict1,f_path+str(cancer)+"_1.txt")
clust.probe_dict_to_txt(dict2,f_path+str(cancer)+"_2.txt")
clust.probe_dict_to_txt(dict3,f_path+str(cancer)+"_3.txt")
clust.probe_dict_to_txt(dict4,f_path+str(cancer)+"_4.txt")


# Generate consensus sequences of length 24,981

sub_cancers = ["brca_1","brca_2","probe_dict_dis","gbm_1","gbm_2","kirc_1","kirc_2","luad_1","luad_2","probe_dict_dis","ov_1","ov_2","ov_3","ov_4","probe_dict_dis"]
cancers = ["brca","brca","coad","gbm","gbm","kirc","kirc","luad","luad","lusc","ov","ov","ov","ov","ucec"]
probes_27_short = open("/Users/rathikannan/Documents/hm_27k_profile/Annotation/probes_27_short.txt","r").readline().strip().split('\t')

f_path = "/Users/rathikannan/Documents/hm_27k_profile/Output/"

# Final matrix containing consensus sequences
# Rows: by sub_cancer
# Cols: by gene index in probes_27
# Will be used to generate hamming distances
all_sub = [None]*len(sub_cancers)

for c in range(0,len(sub_cancers)):
	
	cancer = sub_cancers[c]
	
	c_file = cancers[c]
	
	print "Started: "+str(cancer)
	
	probe_dict = clust.read_probe_dict(f_path+c_file+"/"+cancer+".txt")
	
	final_ls = [None]*24981
	
	# rows index probes
	# columns index samples
	cancer_m = clust.dict_to_matrix_2(probe_dict)
	for i in range(0,len(cancer_m)):
		final_ls[i] = clust.resolve(cancer_m[i])
	
	all_sub[c] = final_ls


	
# Generate hamming distances and write to file
all_sub_hamm = clust.hamming_full_csv(all_sub,f_path+"all_sub_hamming.csv")






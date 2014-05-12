import os
import cPickle as pickle
import pandas as pd
import resolve


hyper_ls = open("hypermeth_targets.txt","r").readline().strip().split('\r')

adf_27 = pickle.load(open("adf_27.p","rb"))

symbol_match = {}
for probe in adf_27.index:
	symb = adf_27.loc[probe,"Symbol"]
	syn = adf_27.loc[probe,"Synonym"]
	if syn == "":
		symbol_match[probe] = [symb]
	else:
		syn = syn.split()
		syn = [i.strip(";") for i in syn]
		syn.append(symb)
		symbol_match[probe] = syn

hyper_dict = {}
for gene in hyper_ls:
	hyper_dict[gene] = []
	
for gene in hyper_ls:
	for probe in symbol_match:
		names = symbol_match[probe]
		if gene in names:
			temp = hyper_dict[gene]
			temp.append(probe)
			hyper_dict[gene] = temp

# Check
print len(hyper_ls)
print len(hyper_dict)
for gene in hyper_dict:
	if hyper_dict[gene] == []:
		print gene
# Genes in hyper_ls that are not found in probes:
#Dc42
#PIZ
#HIN
#COL5A
#p14-ARF

sub_cancers = ["brca_1","brca_2","coad","gbm_1","gbm_2","kirc_1","kirc_2","luad_1","luad_2","lusc","ov_1","ov_2","ov_3","ov_4","ucec"]
cancers = ["brca","coad","gbm","kirc","luad","lusc","ov","ucec"]

temp = hyper_dict.keys()
for gene in temp:
	if len(hyper_dict[gene]) == 0:
		hyper_dict.pop(gene,"None")

# 39 genes represented
hyper_index = sorted(hyper_dict.keys())

probes_27_short = open("probes_27_short.txt").readline().strip().split('\t')

# To build hyper dict per cancer
# Every gene in hyper_dict_i should have associated indices
hyper_dict_i = {}
for gene in hyper_dict:
	indices = []
	for probe in hyper_dict[gene]:
		indices.append(probes_27.index(probe))
	hyper_dict_i[gene] = indices

for cancer in cancers:
	print "started: "+str(cancer)
	cancer_dict = {}
	probe_dict_dis = clust.read_probe_dict(f_path+str(cancer)+"/probe_dict_dis.p")
	for sample in probe_dict_dis:
		betas = probe_dict_dis[sample]
		hypers = []
		for gene in hyper_index:
			indices = hyper_dict_i[gene]
			b = []
			for i in indices:
				b.append(betas[i])
			hypers.append(scored(b))
		cancer_dict[sample] = hypers
	pickle.dump(cancer_dict,open(f_path+"Hyper/"+str(cancer)+".p","wb"))


# To segregate hyper dict into cancer subtypes
cancers = ["brca","gbm","kirc","luad"]
for cancer in cancers:
	
	print "Started: "+str(cancer)
	f_path = "/Users/rathikannan/Documents/Research/Output_27A/"+str(cancer)+"/Cluster/All/"
	f_path2 = "/Users/rathikannan/Documents/Research/Output_27A/Hyper/"
	
	cut1 = open(f_path+"cut1.csv","r").readlines()
	cut1 = [i.strip() for i in cut1]
	cut2 = open(f_path+"cut2.csv","r").readlines()
	cut2 = [i.strip() for i in cut2]

	hyper_dict = pickle.load(open(f_path2+str(cancer)+".p","rb"))
	dict1 = {}
	dict2 = {}
	for key in hyper_dict:
		if key in cut1:
			dict1[key] = hyper_dict[key]
		else:
			dict2[key] = hyper_dict[key]
	
	pickle.dump(dict1,open(f_path2+str(cancer)+"_1.p","wb"))
	pickle.dump(dict2,open(f_path2+str(cancer)+"_2.p","wb"))

cancer = "ov"
print "Started: "+str(cancer)
f_path = "/Users/rathikannan/Documents/Research/Output_27A/"+str(cancer)+"/Cluster/All/"
f_path2 = "/Users/rathikannan/Documents/Research/Output_27A/Hyper/"
	
cut1 = open(f_path+"cut1.csv","r").readlines()
cut1 = [i.strip() for i in cut1]
cut2 = open(f_path+"cut2.csv","r").readlines()
cut2 = [i.strip() for i in cut2]
cut3 = open(f_path+"cut3.csv","r").readlines()
cut3 = [i.strip() for i in cut3]
cut4 = open(f_path+"cut4.csv","r").readlines()
cut4 = [i.strip() for i in cut4]

hyper_dict = pickle.load(open(f_path2+str(cancer)+".p","rb"))
dict1 = {}
dict2 = {}
dict3 = {}
dict4 = {}
for key in hyper_dict:
	if key in cut1:
		dict1[key] = hyper_dict[key]
	elif key in cut2:
		dict2[key] = hyper_dict[key]
	elif key in cut3:
		dict3[key] = hyper_dict[key]
	else:
		dict4[key] = hyper_dict[key]

pickle.dump(dict1,open(f_path2+str(cancer)+"_1.p","wb"))
pickle.dump(dict2,open(f_path2+str(cancer)+"_2.p","wb"))
pickle.dump(dict3,open(f_path2+str(cancer)+"_3.p","wb"))
pickle.dump(dict4,open(f_path2+str(cancer)+"_4.p","wb"))

sub_cancer_hyper_consensus = {}
for cancer in sub_cancers:
	cancer_dict = pickle.load(open(f_path+"Hyper/"+str(cancer)+".p","rb"))
	cancer_df = pd.DataFrame.from_dict(cancer_dict,orient="index")
	cancer_df.columns = hyper_index
	con_ls = []
	for gene in hyper_index:
		temp = np.array(cancer_df[gene])
		con_ls.append(scored(temp))
	sub_cancer_hyper_consensus[cancer] = con_ls

# Write to file
f_out = open(f_path+"Hyper/sub_cancer_hyper_consensus.csv","w")
f_line = ["Sub_Cancers"]+hyper_index
f_out.write("\t".join(f_line)+"\n")
for cancer in sub_cancer_hyper_consensus:
	f_out.write(str(cancer)+"\t"+"\t".join([str(i) for i in sub_cancer_hyper_consensus[cancer]])+"\n")
f_out.close()

# Calculate hamming and write to file
# Cluster by samples
by_samples = dict_to_matrix_1(sub_cancer_hyper_consensus)
hamming_full_csv(by_samples,f_path+"Hyper/sub_cancers_samples_hamming.csv")
by_genes = dict_to_matrix_2(sub_cancer_hyper_consensus)
hamming_full_csv(by_genes,f_path+"Hyper/sub_cancers_genes_hamming.csv")
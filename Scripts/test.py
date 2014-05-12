import os
import cluster_methods as clust

cancer = "brca"

print "Started: "+str(cancer)
	
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
print len(dict1)
print len(dict2)
print len(dict1.values()[0])
print len(dict2.values()[0])

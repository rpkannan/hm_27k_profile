
# Read in distance matrix generated in Python (distance between samples)
sub_dist = as.matrix(read.table("all_sub_hamming.csv",header=FALSE,sep="\t"))
attr(sub_dist,"Size") = nrow(sub_dist)
# Cluster
sub_clust = hclust(sub_dist,method="ward")
plot(sub_clust,labels=rownames(all_sub),main="Association Between Cancer Sub-Types")


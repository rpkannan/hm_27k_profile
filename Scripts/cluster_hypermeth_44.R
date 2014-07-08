setwd("~/Documents/Research/Output_27A/Hyper")
library("gplots")
library(RColorBrewer)
mypal = brewer.pal(3,"RdBu")

# Read data
brca_1 = read.csv("brca_1.csv",header=TRUE,sep="\t")
rownames(brca_1) = brca_1$Samples
brca_1$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
brca_1_dist = as.matrix(read.table("brca_1_hamming_genes.csv",header=FALSE,sep="\t"))
attr(brca_1_dist,"Size") = nrow(brca_1_dist)
#Create heatmap
brca_1_matrix = data.matrix(brca_1)
brca_1_hmap = heatmap.2(brca_1_matrix,Rowv=NA,distfun = function(x) as.dist(brca_1_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="BRCA_1")

# Read data
brca_2 = read.csv("brca_2.csv",header=TRUE,sep="\t")
rownames(brca_2) = brca_2$Samples
brca_2$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
brca_2_dist = as.matrix(read.table("brca_2_hamming_genes.csv",header=FALSE,sep="\t"))
attr(brca_2_dist,"Size") = nrow(brca_2_dist)
#Create heatmap
brca_2_matrix = data.matrix(brca_2)
brca_2_hmap = heatmap.2(brca_2_matrix,Rowv=NA,distfun = function(x) as.dist(brca_2_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="BRCA_2")

# Read data
coad = read.csv("coad.csv",header=TRUE,sep="\t")
rownames(coad) = coad$Samples
coad$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
coad_dist = as.matrix(read.table("coad_hamming_genes.csv",header=FALSE,sep="\t"))
attr(coad_dist,"Size") = nrow(coad_dist)
#Create heatmap
coad_matrix = data.matrix(coad)
coad_hmap = heatmap.2(coad_matrix,Rowv=NA,distfun = function(x) as.dist(coad_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="COAD")


# Read data
gbm_1 = read.csv("gbm_1.csv",header=TRUE,sep="\t")
rownames(gbm_1) = gbm_1$Samples
gbm_1$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
gbm_1_dist = as.matrix(read.table("gbm_1_hamming_genes.csv",header=FALSE,sep="\t"))
attr(gbm_1_dist,"Size") = nrow(gbm_1_dist)
#Create heatmap
gbm_1_matrix = data.matrix(gbm_1)
gbm_1_hmap = heatmap.2(gbm_1_matrix,Rowv=NA,distfun = function(x) as.dist(gbm_1_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="GBM_1")



# Read data
gbm_2 = read.csv("gbm_2.csv",header=TRUE,sep="\t")
rownames(gbm_2) = gbm_2$Samples
gbm_2$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
gbm_2_dist = as.matrix(read.table("gbm_2_hamming_genes.csv",header=FALSE,sep="\t"))
attr(gbm_2_dist,"Size") = nrow(gbm_2_dist)
#Create heatmap
gbm_2_matrix = data.matrix(gbm_2)
gbm_2_hmap = heatmap.2(gbm_2_matrix,Rowv=NA,distfun = function(x) as.dist(gbm_2_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="GBM_2")




# Read data
kirc_1 = read.csv("kirc_1.csv",header=TRUE,sep="\t")
rownames(kirc_1) = kirc_1$Samples
kirc_1$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
kirc_1_dist = as.matrix(read.table("kirc_1_hamming_genes.csv",header=FALSE,sep="\t"))
attr(kirc_1_dist,"Size") = nrow(kirc_1_dist)
#Create heatmap
kirc_1_matrix = data.matrix(kirc_1)
kirc_1_hmap = heatmap.2(kirc_1_matrix,Rowv=NA,distfun = function(x) as.dist(kirc_1_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="KIRC_1")




# Read data
kirc_2 = read.csv("kirc_2.csv",header=TRUE,sep="\t")
rownames(kirc_2) = kirc_2$Samples
kirc_2$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
kirc_2_dist = as.matrix(read.table("kirc_2_hamming_genes.csv",header=FALSE,sep="\t"))
attr(kirc_2_dist,"Size") = nrow(kirc_2_dist)
#Create heatmap
kirc_2_matrix = data.matrix(kirc_2)
kirc_2_hmap = heatmap.2(kirc_2_matrix,Rowv=NA,distfun = function(x) as.dist(kirc_2_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="KIRC_2")





# Read data
luad_1 = read.csv("luad_1.csv",header=TRUE,sep="\t")
rownames(luad_1) = luad_1$Samples
luad_1$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
luad_1_dist = as.matrix(read.table("luad_1_hamming_genes.csv",header=FALSE,sep="\t"))
attr(luad_1_dist,"Size") = nrow(luad_1_dist)
#Create heatmap
luad_1_matrix = data.matrix(luad_1)
luad_1_hmap = heatmap.2(luad_1_matrix,Rowv=NA,distfun = function(x) as.dist(luad_1_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="LUAD_1")




# Read data
luad_2 = read.csv("luad_2.csv",header=TRUE,sep="\t")
rownames(luad_2) = luad_2$Samples
luad_2$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
luad_2_dist = as.matrix(read.table("luad_2_hamming_genes.csv",header=FALSE,sep="\t"))
attr(luad_2_dist,"Size") = nrow(luad_2_dist)
#Create heatmap
luad_2_matrix = data.matrix(luad_2)
luad_2_hmap = heatmap.2(luad_2_matrix,Rowv=NA,distfun = function(x) as.dist(luad_2_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="LUAD_2")




# Read data
lusc = read.csv("lusc.csv",header=TRUE,sep="\t")
rownames(lusc) = lusc$Samples
lusc$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
lusc_dist = as.matrix(read.table("lusc_hamming_genes.csv",header=FALSE,sep="\t"))
attr(lusc_dist,"Size") = nrow(lusc_dist)
#Create heatmap
lusc_matrix = data.matrix(lusc)
lusc_hmap = heatmap.2(lusc_matrix,Rowv=NA,distfun = function(x) as.dist(lusc_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="LUSC")




# Read data
ov_1 = read.csv("ov_1.csv",header=TRUE,sep="\t")
rownames(ov_1) = ov_1$Samples
ov_1$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
ov_1_dist = as.matrix(read.table("ov_1_hamming_genes.csv",header=FALSE,sep="\t"))
attr(ov_1_dist,"Size") = nrow(ov_1_dist)
#Create heatmap
ov_1_matrix = data.matrix(ov_1)
ov_1_hmap = heatmap.2(ov_1_matrix,Rowv=NA,distfun = function(x) as.dist(ov_1_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="OV_1")




# Read data
ov_2 = read.csv("ov_2.csv",header=TRUE,sep="\t")
rownames(ov_2) = ov_2$Samples
ov_2$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
ov_2_dist = as.matrix(read.table("ov_2_hamming_genes.csv",header=FALSE,sep="\t"))
attr(ov_2_dist,"Size") = nrow(ov_2_dist)
#Create heatmap
ov_2_matrix = data.matrix(ov_2)
ov_2_hmap = heatmap.2(ov_2_matrix,Rowv=NA,distfun = function(x) as.dist(ov_2_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="OV_2")




# Read data
ov_3 = read.csv("ov_3.csv",header=TRUE,sep="\t")
rownames(ov_3) = ov_3$Samples
ov_3$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
ov_3_dist = as.matrix(read.table("ov_3_hamming_genes.csv",header=FALSE,sep="\t"))
attr(ov_3_dist,"Size") = nrow(ov_3_dist)
#Create heatmap
ov_3_matrix = data.matrix(ov_3)
ov_3_hmap = heatmap.2(ov_3_matrix,Rowv=NA,distfun = function(x) as.dist(ov_3_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="OV_3")



# Read data
ov_4 = read.csv("ov_4.csv",header=TRUE,sep="\t")
rownames(ov_4) = ov_4$Samples
ov_4$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
ov_4_dist = as.matrix(read.table("ov_4_hamming_genes.csv",header=FALSE,sep="\t"))
attr(ov_4_dist,"Size") = nrow(ov_4_dist)
#Create heatmap
ov_4_matrix = data.matrix(ov_4)
ov_4_hmap = heatmap.2(ov_4_matrix,Rowv=NA,distfun = function(x) as.dist(ov_4_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="OV_4")




# Read data
ucec = read.csv("ucec.csv",header=TRUE,sep="\t")
rownames(ucec) = ucec$Samples
ucec$Samples = NULL
# Read in distance matrix generated in Python (distance between genes)
ucec_dist = as.matrix(read.table("ucec_hamming_genes.csv",header=FALSE,sep="\t"))
attr(ucec_dist,"Size") = nrow(ucec_dist)
#Create heatmap
ucec_matrix = data.matrix(ucec)
ucec_hmap = heatmap.2(ucec_matrix,Rowv=NA,distfun = function(x) as.dist(ucec_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,labRow=NA,xlab="Genes",main="UCEC")

# Read data
all_sub = read.csv("sub_cancer_hyper_consensus.csv",header=TRUE,sep="\t")
rownames(all_sub) = all_sub$Sub_Cancers
all_sub$Sub_Cancers = NULL
# Read in distance matrix generated in Python (distance between samples)
all_sub_dist = as.matrix(read.table("sub_cancers_samples_hamming.csv",header=FALSE,sep="\t"))
attr(all_sub_dist,"Size") = nrow(all_sub_dist)
#Create heatmap
all_sub_matrix = data.matrix(all_sub)
all_sub_hmap = heatmap.2(all_sub_matrix,Colv=NA,distfun = function(x) as.dist(all_sub_dist), hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,xlab="Genes",ylab="Cancer Sub-Types",main="By Cancer Sub-Type")
all_sub_test = heatmap.2(all_sub_matrix, hclustfun = function(x) hclust(x,method = 'ward'),col=mypal,trace="none",scale="none",na.rm = TRUE, key=TRUE,xlab="Genes",ylab="Cancer Sub-Types",main="By Cancer Sub-Type")




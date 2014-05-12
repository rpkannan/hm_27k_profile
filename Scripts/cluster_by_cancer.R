cancers= c("brca","coad","gbm","kirc","luad","lusc","ov","ucec")
f_path = getwd()+"/Output/"

library("gplots")
library(RColorBrewer)
mypal = brewer.pal(3,"RdBu")


# BRCA
# Read in distance matrix generated in Python (distance between samples)
brca_dist = as.matrix(read.table(paste(f_path,"brca/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(brca_dist,"Size") = nrow(brca_dist)
# Cluster
brca_clust = hclust(brca_dist,method="ward")
png(filename=paste(f_path,"brca/by_sample_dendro.png"))
plot(brca_clust,labels=FALSE,main="BRCA (n=315)")
rect.hclust(brca_clust, k=2, border="blue")
dev.off()
brca_cut = cutree(brca_clust,k=2) # 21117.400
brca_samples = data.frame(matrix(NA, nrow = 315, ncol = 2))
brca_samples$X1 = rownames(brca_data)
brca_samples$X2 = brca_cut
brca_cut1 = brca_samples$X1[brca_samples$X2 == 1]
write.table(brca_cut1,file=paste(f_path,"brca/cut1.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
brca_cut2 = brca_samples$X1[brca_samples$X2 == 2]
write.table(brca_cut2,file=paste(f_path,"brca/cut2.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


# COAD
# Read in distance matrix generated in Python (distance between samples)
coad_dist = as.matrix(read.table(paste(f_path,"coad/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(coad_dist,"Size") = nrow(coad_dist)
# Cluster
coad_clust = hclust(coad_dist,method="ward")
png(filename=paste(f_path,"coad/by_sample_dendro.png"))
plot(coad_clust,labels=FALSE,main="COAD (n=166)") # 15563.879
dev.off()


# GBM
# Read in distance matrix generated in Python (distance between samples)
gbm_dist = as.matrix(read.table(paste(f_path,"gbm/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(gbm_dist,"Size") = nrow(gbm_dist)
# Cluster
gbm_clust = hclust(gbm_dist,method="ward")
png(filename=paste(f_path,"gbm/by_sample_dendro.png"))
plot(gbm_clust,labels=FALSE,main="GBM (n=295)")
rect.hclust(gbm_clust, k=2, border="blue") # 21342.033
dev.off()
gbm_cut = cutree(gbm_clust,k=2)
gbm_samples = data.frame(matrix(NA, nrow = 295, ncol = 2))
gbm_samples$X1 = rownames(gbm_data)
gbm_samples$X2 = gbm_cut
gbm_cut1 = gbm_samples$X1[gbm_samples$X2 == 1]
write.table(gbm_cut1,file=paste(f_path,"gbm/cut1.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
gbm_cut2 = gbm_samples$X1[gbm_samples$X2 == 2]
write.table(gbm_cut2,file=paste(f_path,"gbm/cut2.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


# KIRC
# Read in distance matrix generated in Python (distance between samples)
kirc_dist = as.matrix(read.table(paste(f_path,"kirc/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(kirc_dist,"Size") = nrow(kirc_dist)
# Cluster
kirc_clust = hclust(kirc_dist,method="ward")
png(filename=paste(f_path,"kirc/by_sample_dendro.png"))
plot(kirc_clust,labels=FALSE,main="KIRC (n=219)")
rect.hclust(kirc_clust, k=2, border="blue") # 19152.379
dev.off()
kirc_cut = cutree(kirc_clust,k=2)
kirc_samples = data.frame(matrix(NA, nrow = 219, ncol = 2))
kirc_samples$X1 = rownames(kirc_data)
kirc_samples$X2 = kirc_cut
kirc_cut1 = kirc_samples$X1[kirc_samples$X2 == 1]
write.table(kirc_cut1,file=paste(f_path,"kirc/cut1.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
kirc_cut2 = kirc_samples$X1[kirc_samples$X2 == 2]
write.table(kirc_cut2,file=paste(f_path,"kirc/cut2.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


# LUAD
# Read in distance matrix generated in Python (distance between samples)
luad_dist = as.matrix(read.table(paste(f_path,"luad/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(luad_dist,"Size") = nrow(luad_dist)
# Cluster
luad_clust = hclust(luad_dist,method="ward")
png(filename=paste(f_path,"luad/by_sample_dendro.png"))
plot(luad_clust,labels=FALSE,main="LUAD (n=126)")
rect.hclust(luad_clust, k=2, border="blue") # 21050.091
dev.off()
luad_cut = cutree(luad_clust,k=2)
luad_samples = data.frame(matrix(NA, nrow = 126, ncol = 2))
luad_samples$X1 = rownames(luad_data)
luad_samples$X2 = luad_cut
luad_cut1 = luad_samples$X1[luad_samples$X2 == 1]
write.table(luad_cut1,file=paste(f_path,"luad/cut1.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
luad_cut2 = luad_samples$X1[luad_samples$X2 == 2]
write.table(luad_cut2,file=paste(f_path,"luad/cut1.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


# LUSC
# Read in distance matrix generated in Python (distance between samples)
lusc_dist = as.matrix(read.table(paste(f_path,"lusc/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(lusc_dist,"Size") = nrow(lusc_dist)
# Cluster
lusc_clust = hclust(lusc_dist,method="ward")
png(filename=paste(f_path,"lusc/by_sample_dendro.png"))
plot(lusc_clust,labels=FALSE,main="LUSC (n=134)") # 16895.093
dev.off()


# OV
# Read in distance matrix generated in Python (distance between samples)
ov_dist = as.matrix(read.table(paste(f_path,"ov/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(ov_dist,"Size") = nrow(ov_dist)
# Cluster
ov_clust = hclust(ov_dist,method="ward")
png(filename=paste(f_path,"ov/by_sample_dendro.png"))
plot(ov_clust,labels=FALSE,main="OV (n=591)")
rect.hclust(ov_clust, k=4, border="blue") # 19306.387 22028.307 34561.637
dev.off()
ov_cut = cutree(ov_clust,k=4)
ov_samples = data.frame(matrix(NA, nrow = 591, ncol = 2))
ov_samples$X1 = rownames(ov_data)
ov_samples$X2 = ov_cut
ov_cut1 = ov_samples$X1[ov_samples$X2 == 1]
write.table(ov_cut1,file=paste(f_path,"ov/cut1.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
ov_cut2 = ov_samples$X1[ov_samples$X2 == 2]
write.table(ov_cut2,file=paste(f_path,"ov/cut2.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
ov_cut3 = ov_samples$X1[ov_samples$X2 == 3]
write.table(ov_cut3,file=paste(f_path,"ov/cut3.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
ov_cut4 = ov_samples$X1[ov_samples$X2 == 4]
write.table(ov_cut4,file=paste(f_path,"ov/cut4.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)


# UCEC
# Read in distance matrix generated in Python (distance between samples)
ucec_dist = as.matrix(read.table(paste(f_path,"ucec/probe_dis_hamming_samples.txt",sep=""),header=FALSE,sep="\t"))
attr(ucec_dist,"Size") = nrow(ucec_dist)
# Cluster
ucec_clust = hclust(ucec_dist,method="ward")
png(filename=paste(f_path,"ucec/by_sample_dendro.png"))
plot(ucec_clust,labels=FALSE,main="UCEC (n=117)") #13670.358
dev.off()


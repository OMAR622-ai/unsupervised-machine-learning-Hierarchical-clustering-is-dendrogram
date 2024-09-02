#Load libraries
library(stats)
library(tools)

#Load data
data <- read.table ("https://raw.githubusercontent.com/pine-bio-support/Unsupervised-codeomics/main/CellLines_7000Genes_marked_UI.txt", header = TRUE, row.names = 1)
#Check dimensions of data
dim(data)

#Transpose the data
data1 <- t(data)

#Check dimensions of data
dim(data1)

# Define clustering and save in a object 
#specify distance matrix (euclidean / manhattan / maximum / canberra / binary / minkowski) and linkage(single / complete / average / mean / centroid / ward.D / ward.D2) type
cl <- hclust(dist (data1, method="euclidean"), method="complete")

#Plot clusters
plot(cl, cex = 0.6)

#Define parameters and save in object
descr1 <- paste ("Distance: ", "euclidean", sep="")
descr2 <- paste ("Linkage: ", "complete", sep="")

#Plot Clusters
plot(cl, xlab=descr1, sub=descr2, cex = 0.6)

#cutree to get the required number of clusters and output cluster IDs for all the samples
clust_ids <- cutree(cl, k=5)
head(clust_ids)

#Output table with cluster IDs 
txt_file_c <- paste("file_name", "_", "cRes_hclust.ClustersOnly", ".txt", sep="")
cat ("Object\tClusterID\n", file=txt_file_c, append=FALSE)


#Write into a file
write.table (clust_ids, file=txt_file_c, append=TRUE, sep="\t", col.names=FALSE, row.names=TRUE, quote = FALSE)


#Output full table with cluster IDs 
txt_file <- paste("file_name", "_", "cRes_hclust.FullTable", ".txt", sep="")	
out_table <- cbind (clust_ids, data)
colnames (out_table) [1] <- "ClusterID"

#Write into a file
write.table (out_table, file=txt_file, append=TRUE, sep="\t", col.names=TRUE, row.names=TRUE)
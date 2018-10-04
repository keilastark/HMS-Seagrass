library(tidyverse)
library(vegan)
library(viridis)
library(gplots)

com_data<-read.csv("./Data/for_ems.csv")

head(com_data)
heatmap.2(as.matrix(vegdist(decostand(com_data[,-1], method = "hellinger"), method = "euclidean", diag = TRUE, upper = TRUE)),col = viridis(100),linecol = NA, trace = "none",density.info = "none")

rownames(com_data) <- com_data[,1]
com_data <- com_data[,-1]

com_clusters<-hclust(vegdist(decostand(com_data, method = "hellinger"),method = "euclidean"),method = "ward.D")
plot(com_clusters)
groups <- cutree(com_clusters, k=3)
rect.hclust(com_clusters, k=3, border="red")


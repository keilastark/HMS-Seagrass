library(tidyverse)
library(vegan)
library(viridis)
library(gplots)

com_data<-read.csv("./Data/for_ems.csv")
blues<-read.csv("./Data/blues_only2.csv")[,-1]
blues<-colnames(blues[,colSums(blues)>0])

purples<-read.csv("./Data/purples_only2.csv")[,-1]
purples<-colnames(purples[,colSums(purples)>0])

greys<-read.csv("./Data/grays_only2.csv")[,-1]
greys<-colnames(greys[,colSums(greys)>0])

species_type <- data.frame(species = c(blues, purples, greys), type = c(rep("blue",length(blues)), rep("purple", length(purples)), rep("grey", length(greys))))

head(com_data)
heatmap.2(as.matrix(vegdist(decostand(com_data[,-1], method = "hellinger"), method = "euclidean", diag = TRUE, upper = TRUE)),col = viridis(100),linecol = NA, trace = "none",density.info = "none")

rownames(com_data) <- com_data[,1]
com_data <- com_data[,-1]

com_clusters<-hclust(vegdist(decostand(com_data, method = "hellinger"),method = "euclidean"),method = "ward.D")
plot(com_clusters)
groups <- cutree(com_clusters, k=3)
rect.hclust(com_clusters, k=3, border="red")

com_stand<-decostand(com_data, margin = 1, method = "total")
plot(apply(com_stand,2,median)[order(apply(com_stand,2,median),decreasing = T)])

species_type_sub <- left_join(data.frame(species=colnames(as.matrix(com_stand)[,order(apply(com_stand,2,mean),decreasing = T)][,1:10])),
                              species_type)
library(heatmaply)
options(scipen=1)
sub_com_ordered<-round(as.matrix(com_stand)[,order(apply(com_stand,2,mean),decreasing = T)][,1:10],digits = 2)

heatmaply(sub_com_ordered, k_row = 3, k_col = 3,Colv = F,hclust_method = "ward.D",)

sub_com_ordered2<-decostand(com_data, method = "hellinger")

heatmaply(sub_com_ordered2, k_row = 3, k_col = 5,hclust_method = "ward.D")



library(Hmsc)
library(tidyverse)

Y <- read.csv("./Data/Y_matrix.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(Y) <- Y[,1]
X <- read.csv("./Data/X_matrix.csv", header = TRUE, stringsAsFactors = FALSE) # Environmental covariate matrix (X)
#X <- X[,-c(7,8,10)] #watch out! this is way of removing data can easily lead to mistakes if your column order changes
X <- X %>% dplyr::select(-dissox) # this is a much safer way to remove data, prevents things messing up if the order of the columns changes
rownames(X) <- Y[,1]

spatial <- read.csv("./Data/coords.csv", header = TRUE, stringsAsFactors = FALSE) # lat&long coordinates
rownames(spatial) <- Y[,1]
#spatial <- spatial[,-1] #why are you removing longitude, this would assume no spatial distance between sites at the same latitude but at different longitudes

studyDesign <- read.csv("./Data/pi.csv", header = TRUE, stringsAsFactors = FALSE) #factors

Y <- Y[,-1] #remove actual site names from Y matrix now that rownames for all datasets have been established

## Remove NA's / get into right format
X[is.na(X)] <- 0
Y[is.na(Y)] <- 0
Y <- as.matrix(Y)
X <- as.data.frame(X)
spatial <- as.data.frame(spatial)
spatial$longitude <- spatial$longitude-min(spatial$longitude)

#spat <- data.frame(spat = sprintf('spatial_%.2d',1:78)) #spatial factor column for studyDesign
#studyDesign <- cbind(studyDesign, spat)
studyDesign$Spatial <- factor(studyDesign$Quadrat)

rL1 = HmscRandomLevel(units = studyDesign$Quadrat)
rL2 = HmscRandomLevel(units = studyDesign$Site)
rL3 = HmscRandomLevel(units = studyDesign$Region)
rL4= HmscRandomLevel(sData = spatial)

# Construct and fit HMSC model
#hM <- Hmsc(Y = Y, XData = X, 
           XFormula = ~eelgrass_lai + ph + sstmean + nitrate + chlomean + salinity, 
           studyDesign = studyDesign, 
           ranLevels = list(Quadrat = rL1, Site = rL2, Region = rL3, Spatial = rL4)) # I have selected only these variables because some covary (ggpairs analysis)

hM <- Hmsc(Y = Y, XData = X, 
           XFormula = ~eelgrass_lai + ph + sstmean + nitrate + chlomean + salinity, 
           studyDesign = studyDesign, 
           ranLevels = list(Spatial = rL4))

mod <- sampleMcmc(hM, samples = 1000 , transient = 1000, thin = 100, verbose = 20000)



preds <- computePredictedValues(mod)
MF = evaluateModelFit(hM=mod, predY=preds)
hist(MF$R2,xlim = c(0,1),main=paste0("Mean = ",round(mean(MF$R2),2)))

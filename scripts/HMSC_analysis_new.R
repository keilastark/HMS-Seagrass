### New HMSC Analysis
library(devtools)
install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))
library(Hmsc)
setwd("~/Documents/seagrass_metacom_paper/FINAL SUBMISSION")

Y <- read.csv("Y_matrix.csv", header = TRUE, stringsAsFactors = FALSE) # Site-by-sp matrix (Y)
rownames(Y) <- Y[,1]

X <- read.csv("X_matrix.csv", header = TRUE, stringsAsFactors = FALSE) # Environmental covariate matrix (X)
rownames(X) <- Y[,1]

spatial <- read.csv("coords.csv", header = TRUE, stringsAsFactors = FALSE) # lat&long coordinates
rownames(spatial) <- Y[,1]
spatial <- spatial[,2:3]

Pi <- read.csv("pi.csv", header = TRUE, stringsAsFactors = FALSE)


Y <- as.matrix(Y)
X <- as.data.frame(X)

Ysub <- Y[c(37,49,55,61),]
Xsub <- X[c(37,49,55,61),]
studyDesignsub <- studyDesign[c(37,49,55,61),]
spatialsub <- spatial[c(37,49,55,61),]

studyDesign <- Pi
spat <- data.frame(spat = sprintf('spatial_%.2d',1:78))
studyDesign <- cbind(studyDesign, spat)

rL1 = HmscRandomLevel(units = studyDesign$Quadrat)
rL2 = HmscRandomLevel(units = studyDesign$Site)
rL3 = HmscRandomLevel(units = studyDesign$Region)
rL4 = HmscRandomLevel(sData = spatial[,2:3])

rL1sub = HmscRandomLevel(units = studyDesignsub$Quadrat)
rL2sub = HmscRandomLevel(units = studyDesignsub$Site)
rL3sub = HmscRandomLevel(units = studyDesignsub$Region)
rL4sub = HmscRandomLevel(sData = spatialsub[,2:3])

library(sp)
library(rgdal)

xy <- data.frame(ID = 1:78, X = spatial[,1], Y = spatial[,2])
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example

res <- spTransform(xy, CRS("+proj=utm +zone=51 ellps=WGS84"))
res

msub = Hmsc(Y = Ysub, XData = Xsub, XFormula = ~eelgrass_lai + detritus_afdm + algae_biomass + fetch + dissox + salinity + ph + sstmean + nitrate + chlomean, studyDesign = studyDesignsub, ranLevels = list(Quadrat = rL1sub, Site = rL2sub, Region = rL3sub))

hM <- Hmsc(Y = Y, XData = X, XFormula = ~eelgrass_lai + dissox + salinity + ph + sstmean + nitrate + chlomean, studyDesign = studyDesign, ranLevels = list(Quadrat = rL1, Site = rL2, Region = rL3, spat = rL4))

mod <- sampleMcmc(hM, samples = 100 , transient = 1000, thin = 100)



for(r in seq_len(hM$nr)){
  if(hM$rL[[r]]$sDim > 0){
    alphapw = hM$rL[[r]]$alphapw
    np = hM$np[r]
    alphaN = nrow(alphapw)
    }
    if(is.null(hM$rL[[r]]$distMat)){
      s = hM$rL[[4]]$s[levels(hM$dfPi[,4]),]
      distance = as.matrix(dist(s))
    } else{
      distance = hM$rL[[4]]$distMat[levels(hM$dfPi[,4]),levels(hM$dfPi[,4])]
    }
    Wg = array(NA, c(np,np,alphaN))
    iWg = array(NA, c(np,np,alphaN))
    RiWg = array(NA, c(np,np,alphaN))
    detWg = rep(NA, alphaN)
    for(ag in 1:alphaN){
      alpha = alphapw[ag,1]
      if(alpha==0){
        W = diag(np)
      } else{
        W = exp(-distance/alpha)
      }
      RW = chol(W)
      iW = chol2inv(RW)
      Wg[,,ag] = W
      iWg[,,ag] = iW
      RiWg[,,ag] = chol(iW)
      detWg[ag] = 2*sum(log(diag(RW)))
    }


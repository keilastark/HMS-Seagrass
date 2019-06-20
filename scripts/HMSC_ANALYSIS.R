### Re-analyzing BC seagrass data with HMSC
### Keila Stark, 2018
### Full analysis: includes figures that were not included in manuscript

setwd("~/Documents/seagrass_metacom_paper/FINAL SUBMISSION")



#### HMSC ANALYSIS ------------------------------------------------------------

### 1. INSTALL AND LOAD REQUIRED PACKAGES
install.packages('devtools')
install.packages('Rcpp')
install.packages('RcppArmadillo')
install.packages('coda')
install.packages('beanplot')
install.packages('circlize')
install.packages('corrplot')
install.packages('coda')
install_github('guiblanchet/HMSC')

library(devtools)
library(HMSC)
library(vegan)
library(tidyverse)
library(viridis)


### 2. READ IN AND PREPARE DATA

Y_full <- read.csv("Y_matrix_full.csv", header = TRUE, stringsAsFactors = FALSE) # Site-by-sp matrix (Y)
rownames(Y) <- Y[,1]

Y <- read.csv("Y_matrix.csv", header = TRUE, stringsAsFactors = FALSE) # Site-by-sp matrix (Y)
rownames(Y) <- Y[,1]

X <- read.csv("X_matrix.csv", header = TRUE, stringsAsFactors = FALSE) # Environmental covariate matrix (X)
rownames(X) <- Y[,1]

spatial <- read.csv("coords.csv", header = TRUE, stringsAsFactors = FALSE) # lat&long coordinates
rownames(spatial) <- Y[,1]

Pi <- read.csv("pi.csv", header = TRUE, stringsAsFactors = FALSE) #factors representing random effects
Pi <- data.frame(apply(Pi,2,as.factor))

Y <- Y[,-1] #remove actual site names from Y matrix

## Remove NA's get into right format for as.HMSCdata to work
X[is.na(X)] <- 0
Y[is.na(Y)] <- 0
Y <- as.matrix(Y)
X <- as.matrix(X)
spatial <- as.data.frame(spatial)
spatial <- cbind(Pi[,1], spatial)




#### 3. CREATE HMSCdata OBJECT WITH THE ABOVE 4 INPUT MATRICES, GENERATE MODEL --------
#Create HMSCdata object with the above four input matrices
hmsc1 <- as.HMSCdata(Y = Y, X = X, Random = Pi, Auto = spatial, interceptX = FALSE)
envonly <- as.HMSCdata(Y = Y, X = X, Random = Pi, interceptX = FALSE)
spaceonly <- as.HMSCdata(Y = Y, Random = Pi, Auto = spatial, interceptX = FALSE)

# set prior distribution
hmscprior <- as.HMSCprior(hmsc1)
hmscparam <- as.HMSCparam(hmsc1, hmscprior)
envprior <- as.HMSCprior(envonly)
envparam <- as.HMSCparam(envonly,envprior)
spaceprior <- as.HMSCprior(spaceonly)
spaceparam <- as.HMSCparam(spaceonly, spaceprior)


## Generate the hierarchical joint species distribution model with MCMC. We used 200000 interations, as the developers of the method note that Poisson distributions sometimes require 10-100 times more iterations to accurately estimate the model's parameters
model <- hmsc(hmsc1, family = "poisson", niter = 200000, nburn = 100000,
              thin = 100)
envmod <-  hmsc(envonly, family = "poisson", niter = 200000, nburn = 100000,
                thin = 100)
spacemod <-  hmsc(spaceonly, family = "poisson", niter = 200000, nburn = 100000,
                  thin = 100)


### 3. Produce mixing, trace, density plots to look at model performance
### Note: None of the figures/ results from this section are included in the paper.

##Producing MCMC traceplots
mixing <- as.mcmc(model, parameters = "paramX") #Mixing object

### Draw trace and density plots for all combination of parameters
plot(mixing)

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing)

### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### True values
truth <- as.vector(hmsc1$param$paramX)
### Average
average <- apply(model$results$estimation$paramX, 1:2, mean)
### 95% confidence intervals
CI.025 <- apply(model$results$estimation$paramX, 1:2, quantile, probs=0.025)
CI.975 <- apply(model$results$estimation$paramX, 1:2, quantile, probs=0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

### Draw confidence interval plots
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", ylab = "", main="paramX")
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch = 19)
### Summary table
paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")
write.csv(paramXCITable, file='beta.csv')



#### 4. COMPUTE MODEL FIT (R2)
### Use Duarte's code
# Poisson deviance
# Arguments:
# Y - the observations
# mu - the predicted means (of Poisson distribution)
Poisson.deviance <- function(Y, mu)
{
  2*sum(ifelse(Y == 0, 
               -(Y-mu), # when the observation is 0
               Y*log(Y/mu) - (Y-mu))) # else
}

PD <- Poisson.deviance(Y, pred.hmsc)



# Explained deviance calculated for each species, then averaged


D2.Poisson <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  mu <- as.matrix(Y.pred) # we assume that the predicted values are the Poisson mu
  
  resid.devs <- null.devs <- numeric(ncol(Y))
  for(i in 1:ncol(Y))
  {
    resid.devs[i] <- Poisson.deviance(Y[,i], mu[,i])
    null.devs[i]  <- Poisson.deviance(Y[,i], mean(Y[,i]))
  }
  
  D2s <- 1 - resid.devs/null.devs
  # change negative values to 0
  
  
  return(D2s)
}


# Use with HMSC
# Y is the species matrix
# pred.hmsc is the predicted species matrix
pred.hmsc <- predict(model)
pred.hmsc <- read.csv("predsite.csv", header = TRUE, stringsAsFactors = FALSE)
D2 <- D2.Poisson(Y, pred.hmsc)
D2.avg <- mean(D2) # average community level D2 value
D2 <- matrix(,nrow = 1, ncol = 41)
colnames(D2) <- colnames(Y)
quad.scale.D2 <- D2

#### 4. VARIATION PARTITIONING OF SPACE, RANDOM EFFECTS AND ENVIRONMENTAL VARIABLES (FIG 3)----------

# First conduct variation partitioning/ create object summarizing result for global model...
varpart <- variPart(model,  c(rep("Shoot density",1),rep("Leaf area index",1),rep("Seagrass AFDM",1), rep("Detritus AFDM",1),rep("Algal AFDM",1),  rep("Fetch",1),rep("Dissolved O2",1),rep("Salinity",1),rep("pH",1),rep("Mean Temp",1), rep("Nitrates",1),rep("Chl a",1)), type = "I")

# Environment only
varpartenv <- variPart(envmod, c(rep("Habitat structure",2),rep("Food availability",3), rep("Fetch",1),rep("Dissolved O2",1),rep("Salinity",1),rep("pH",1),rep("Mean Temp",1), rep("Nitrates",1),rep("Phosphates",1), rep("Silicates",1),rep("Chl a",1)))

#Space only
varpartspatial<- variPart(spacemod, c(rep("intercept",1)))

spatmeans <- colMeans(varpartspatial) #take the mean variance explained by i.e. spatial variables for all species
envmeans <- colMeans(varpartenv) 
wholemeans <- colMeans(varpart)
## Use these values to calculate different fractions of variation explained 
(See "manual_varpart.xlsx" for details)

#Now we are going to create Fig. 3- plotting the Global Model only
# Set interesting graphical parameters to fit all of our species labels below and legend to the right
par(mar=c(12,4,2,12.6))
par(xpd = FALSE)

# Create the figure. The spacing scheme is to separate broad taxonomic groupings.
colnames(varpart) <- c("Shoot density","Leaf area Index", "Eelgrass biomass", "Macroalgae biomass","Detritus biomass","Fetch","Dissolved O2","Salinity","pH","Avg Temp","Nitrates","Chl a","Quadrat","Site","Region","Spatial distance")  

spnames <- rownames(varpart)

spl <- strsplit(sub("_", ".", spnames), ".", fixed = TRUE)
sapply(spl, "[", 2)

gsub("^[^_]*_|\\.[^.]*$", "", spnames)

3) gsubfn::strapplyc extract everything between underscore and dot.

library(gsubfn)
strapplyc(x, "_(.*)\\.", simplify = TRUE)

barplot(t(varpart), las=2, cex.names=0.8, cex.axis=0.75,  col = c("yellow2", "yellow4","springgreen2","springgreen3", "springgreen4" ,"#31688EE6", "#31688ED9", "#31688ECC", "#31688EBF", "#31688EB3", "#31688EA6", "#31688E99", "#440154FF", "#44015480", "#440154B3", "black"), legend.text=paste(colnames(varpart)," ", "(",
                          signif(100*colMeans(varpart),2),"%)", sep=""),
        args.legend=list(yjust = 0.91,xjust=0.01,horiz=F, bty="n",cex=0.75, pt.cex = 0.5, y.intersp = 0.9, x.intersp = 0.4), space = c(0,0,0,0,0,0,0,0,0,0.6,0,0,0,0.6,0,0.6,0,0,0,0,0,0,0,0,0.6,0,0,0.6,0.6,0.6,0.6,0,0,0,0.6,0.6,0,0,0.6,0,0.6), ylab = "Proportion of variance")


#### 5. PLOT DISTANCE-DECAY OF COMMUNITY SIMILARITY (FIG 2) ---------------------
npred <- 100
predAll <- array(NA,dim=list(nrow(hmsc1$Y),ncol(hmsc1$Y),npred))
predAll2 <- array(NA,dim=list(nrow(hmsc1$Y),ncol(hmsc1$Y),npred))
for (i in 1:npred) {
  predAll[,,i] <- predict(model, newdata=hmsc1)
  predAll2[,,i] <- predict(model, newdata=hmsc11)
}
pred <- apply(predAll,c(1,2),mean)
pred2 <- apply(predAll2,c(1,2),mean)
model <- model3
model$data$X <- hmsc1$X
model$data$Y <- hmsc1$Y
# Calculate and plot similarity
similar <- similarity(model3)
similar2 <- kmeans(similar, 13)
XYall <- hmsc1$Auto[[1]][,2:3]
colnames(XYall) <- c("Longitude", "Latitude")
colrs <- brewer.pal(13, "Accent")
plot(XYall,type="n", xlim = c(-132,-123), ylim=c(48, 53))
for (i in unique(similar2$cluster)) {
  points(XYall[which(similar2$cluster==i),],col=colrs[i],pch=15,cex=1)
}


#### 5. SPECIES-TO-SPECIES CO-OCCURRENCE MATRIX (FIG 3)-----------------------------------------

### Plot random effect estimation through correlation matrix
corMat <- corRandomEff(model,cor=TRUE)
dimnames(corMat)

### Isolate the values of interest
ltri <- lower.tri(apply(corMat[,,,1], 1:2, quantile, probs=0.025),
                  diag=TRUE)
### True values
truth <- as.vector(tcrossprod(simulParamEx1$param$paramLatent[[2]])[ltri])
### Average
average <- as.vector(apply(corMat[, , , 2], 1:2, mean)[ltri])
### 95% confidence intervals
corMat.025 <- as.vector(apply(corMat[, , , 2], 1:2, quantile,
                              probs = 0.025)[ltri])
corMat.975 <- as.vector(apply(corMat[, , , 2], 1:2, quantile,
                              probs = 0.975)[ltri])
CI <- cbind(corMat.025, corMat.975)
### Plot the results
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI, truth), type = "n",
     xlab = "", main = "cov(paramLatent[[1,2]])")
abline(h = 0, col = "grey")
arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)
points(1:nrow(CI), truth, col = "red", pch = 19)
### Mixing object
mixing <- as.mcmc(model, parameters = "paramLatent")
### Draw trace and density plots for all combination of paramters
plot(mixing[[2]])
### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing[[2]])
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)

averageCor1 <- averageCor

averageCor1 <- averageCor1[c(2,7,13,14,25,26,27,29,30,35,37,38,42,43,44,46,48,49,50,54),c(2,7,13,14,25,26,27,29,30,35,37,38,42,43,44,46,48,49,50,54)] #2Take 0 most abundant



#### Plot Fig. 4a 
library(corrplot)
corMat <- corRandomEff(model, cor =TRUE)
averageCor <- apply(corMat$Random[ ,,,1], 1:2, mean)
averageCor1 <- averageCor
averageCor1[-0.5 < averageCor1 & averageCor1< 0.5] <- 0

averagesite<- apply(corMat$Random[ ,,,2], 1:2, mean)
averagesite1 <- averagesite
averagesite1[-0.5 < averagesite1 & averagesite1< 0.5] <- 0

averageregion <- apply(corMat$Random[ ,,,3], 1:2, mean)
averageregion1 <- averageregion
averageregion1[-0.5 < averageregion1 & averageregion1< 0.5] <- 0

par(mfrow = c(1,1))

averageCor1 <- averageCor[c(2,7,13,14,25,26,27,29,30,35,37,38,42,43,44,46,48,49,50,54),c(2,7,13,14,25,26,27,29,30,35,37,38,42,43,44,46,48,49,50,54)] # Take 20 most abundant sp
corrplot(averageCor1, method = "color", col = colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method = "ward.D2", tl.col = "black", tl.cex = 0.8)
corrplot(averageregion1, method = "color", col = colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method = "ward.D2", tl.col = "black", tl.cex = 0.8)
corrplot(averagesite1, method = "color", col = colorRampPalette(c("#21908C","white","#440154"))(200), order = "hclust", hclust.method = "ward.D2", tl.col = "black", tl.cex = 0.8)



library(RColorBrewer) ### Plot Fig 1b
n <- nrow(hmsc1$X)
hmsc11 <- hmsc1
hmsc11$X <- matrix(rep(colMeans(hmsc11$X),each=n),
                      nrow=n,ncol=ncol(hmsc11$X))

# Make predictions
npred <- 3003

predAll <- array(NA,dim=list(nrow(hmsc1$Y),ncol(hmsc1$Y),npred))
predAll2 <- array(NA,dim=list(nrow(hmsc11$Y),ncol(hmsc11$Y),npred))
for (i in 1:npred) {
  predAll[,,i] <- predict(model, newdata=hmsc1)
  predAll2[,,i] <- predict(model, newdata=hmsc11)
}
pred <- apply(predAll,c(1,2),mean)
pred2 <- apply(predAll2,c(1,2),mean)


# Plot community similarity as a function of distance, based on pairs of different points
nrepls <- 150
ta <- matrix(0,nrow=nrepls,ncol=3)
for (i in 1:nrepls) {
  s1s2 <- sample(1:nrow(hmsc1$X),2)
  s1 <- s1s2[1]
  s2 <- s1s2[2]
  ta[i,1] <- sqrt(sum((hmsc1$Auto[[1]][s1,2:3]-hmsc1$Auto[[1]][s2,2:3])^2))
  ta[i,2] <- cor(pred[s1,], pred[s2,])
  ta[i,3] <- cor(pred[s1,], pred2[s2,])
}
maxd <- max(ta[,1],na.rm=T)

plot.new()
par(mar = c(4,4,1,1))
par(oma = c(4, 1, 1, 1))
wop <- ta[,1]*111
plot(x=wop,y=ta[,2],pch=1,col="firebrick",xlab="Pairwise site distance (km)",ylab="Predicted community similarity")
lines(y=predict(lm(ta[,2]~ta[,1])),x=wop,col="firebrick",  lwd = 3)
points(x=wop,y=ta[,3],pch=2,col="cyan3",xaxt='n',yaxt="n")
lines(y=predict(lm(ta[,3]~ta[,1])),x=wop,col="cyan3",lwd = 3)
legend("bottomright", legend=c("Space + environmental covariates","Space only"),col=c("firebrick","cyan3") ,cex= 0.75, pch=c(1,2), xpd = NA, inset = c(.5,-.5), bty = "n")



plot.new()

# EMS ordination image ----------------------------------------------------
library(metacom)
ems <- read.csv("for_ems.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(ems) <- c("DC", "RB", "SA", "DK", "EB", "IN", "CB", "GB", "JB", "LH", "SS", "HL", "RA")
ems <- ems[,-1]

#Get output of EMS 
ordermatrix_abundance <- OrderMatrix(ems, binary = FALSE)

Imagine(ems, binary = FALSE)

write.csv(ordermatrix_abundance, file = "emsresult_abundance.csv")

emsresult <- read.csv("emsresult_edited.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(emsresult) <- emsresult[,1]
emsresult1 <- emsresult[,2:58]
emsresult1 <- as.numeric(emsresult1)
emsresult1 <- as.data.frame(emsresult1)
emsresult1 <- as.matrix(emsresult1)

emsresult2 <- read.csv("emsresult_abundance_edited.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(emsresult2) <- emsresult2[,1]
emsresult2 <- emsresult2[,2:58]
emsresult2 <- as.numeric(emsresult1)
emsresult2 <- as.data.frame(emsresult1)
emsresult2 <- as.matrix(emsresult1)



### Made a modified version of the Imagine function from the metacom package to fill cells with colour based on which co-occurrence grouping they belonged to (from Ward D2 cluster analysis on HMSC model)


imagineog <- function (comm, col = c("white","black"), fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = col , xlab = "", ylab = "", axes = FALSE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}

imagine2(emsresult1)
imagine2(emsresult2)


gastropods <- Y[,1:12]
bivalves <- Y[,13:17]
polychaetes <- Y[,18:25]
gammarid_amphipods <- Y[,26:41]
caprellid_amphipods <- Y[,42:44]
cumaceans <- Y[,45]
tanaids <- Y[,46]
leptostracans <- Y[,47]
copepods <- Y[,48:51]
pycnogonids <- Y[,52]
isopods <- Y[,53:55]
crabs <- Y[,56:57]
percarids <- Y[,58]

gastro.total <- rowSums(gastropods)
bivalves.total <- rowSums(bivalves)
polychaetes.total <- rowSums(polychaetes)
gammarids.total <- rowSums(gammarid_amphipods)
caprellids.total <- rowSums(caprellid_amphipods)
copepods.total <- rowSums(copepods)
iso.total <- rowSums(isopods)
crab.total <- rowSums(crabs)


stacked <- cbind( gastro.total , bivalves.total , polychaetes.total , gammarids.total ,caprellids.total , cumaceans, tanaids, leptostracans , copepods.total , pycnogonids , iso.total , crab.total , percarids)

stacked_collapsed <- data.frame(HL= colSums(stacked[31:36,]), Ramsay = colSums(stacked[55:60,]), Ducking = colSums(stacked[13:18,]), Elbow = colSums(stacked[19:24,]),  Indian = colSums(stacked[37:42,]), Sarita = colSums(stacked[67:72,]), Robbers= colSums(stacked[61:66,]), Dodger = colSums(stacked[7:12,]),  Sidney = colSums(stacked[73:78,]), Cabbage = colSums(stacked[1:6,]),  Gallagher = colSums(stacked[25:30,]),  James = colSums(stacked[43:48,]), Lyall = colSums(stacked[49:54,]))

stacked3 <- t(stacked_collapsed)

for(i in 1:nrow(stacked3)){
  stacked3[i,] <- stacked3[i,]/sum(stacked3[i,])
}
rowSums(stacked3) #Checking to make sure for-loop made rowsums all equal to 1 for barplot


par(mar = c(5,5,4,10))

par(xpd = TRUE)
barplot(t(stacked3), col = viridis(13), names.arg = rownames(stacked3), space = 0, cex.names = "0.9", xlab = "Site", ylab = "Relative abundance", legend.text = c("gastropods","bivalves","polychaetes","gammarids", "caprellids", "cumaceans" , "tanaids","leptostracans", "copepods","pycnogonids","isopods","crabs","percarids"),  args.legend=list(x="topright", inset = c(-0.47,0), xpd = "TRUE", cex = 0.9, bty = "n", x.intersp = 0.2, title = "Taxonomic group", title.adj = 0.15))

### Trying Patrick's last suggestion
emstransformed <- ems ## create new dataframe and for-loop it to divide species count per site by the maximum abundance. Just to scale so super abundance species don't overwhelm the figure

for(i in 1:nrow(emstransformed)){
  emstransformed[i,] <- emstransformed[i,]/max(emstransformed[i,])
}

emstransformedt <- t(emstransformed)

emstransformed1 <- data.frame( HL = emstransformedt[,12], IN = emstransformedt[,6], RB = emstransformedt[,3], CB = emstransformedt[,7], DK = emstransformedt[,5], DC = emstransformedt[,1], RA = emstransformedt[,13], GB = emstransformedt[,8], EB = emstransformedt[,5], JB = emstransformedt[,9], LH = emstransformedt[,10], SA = emstransformedt[,3], SS = emstransformedt[,11])

write.csv(emstransformed, file = "emstransformed.csv") ### put this into excel to separate out the three corrplot groupings into 3 separate data frames

logtransformed <- ems

for(i in 1:ncol(logtransformed)){
  logtransformed[,i] <- log(logtransformed[,i])
}

logtransformed[logtransformed == -Inf] <- 0

write.csv(logtransformed, file = "logtransformed.csv")

## read them back in after editing in excel 
grays_only3 <- read.csv("grays_only3.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(grays_only3) <- grays_only3[,1]
grays_only3 <- grays_only3[,2:56]


blues_only3 <- read.csv("blues_only3.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(blues_only3) <- blues_only3[,1]
blues_only3 <- blues_only3[,2:56]


purples_only3 <- read.csv("purples_only3.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(purples_only3) <- purples_only3[,1]
purples_only3 <- purples_only3[,2:56]

## Do this again for log transformed stuff 

log_grays_only <- read.csv("log_grays_only.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(log_grays_only) <- log_grays_only[,1]
log_grays_only1 <- log_grays_only[,2:56]


log_blues_only <- read.csv("log_blues_only.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(log_blues_only) <- log_blues_only[,1]
log_blues_only1 <- log_blues_only[,2:56]


log_purples_only <- read.csv("log_purples_only.csv", header = TRUE, stringsAsFactors = FALSE)
row.names(log_purples_only) <- log_purples_only[,1]
log_purples_only1 <- log_purples_only[,2:56]






 #### black and white Imagein function
imagine1 <- function (comm, col = c(0,1), fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = c("white", "grey50"), xlab = "", ylab = "", axes = FALSE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}


ns <- read.csv("n_s_matrix.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(ns) <- ns[,1]
ns <- ns[,-1]
ns[ns > 1] <- 1
imagine1(ns)



## make new imagine functions for each of the data frames
imagine3 <- function (comm, col = col, fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = c("#00000000","#52525203","#52525205","#52525208","#5252520A","#5252520D","#5252520F","#52525212","#52525214","#52525217","#5252521A","#5252521C","#5252521F","#52525221","#52525224","#52525226","#52525229","#5252522B","#5252522E","#52525230","#52525233","#52525236","#52525238","#5252523B","#5252523D","#52525240","#52525242","#52525245","#52525247","#5252524A","#5252524D","#5252524F","#52525252","#52525254","#52525257","#52525259","#5252525C","#5252525E","#52525261","#52525263","#52525266","#52525269","#5252526B","#5252526E","#52525270","#52525273","#52525275","#52525278","#5252527A","#5252527D","#52525280","#52525282","#52525285","#52525287","#5252528A","#5252528C","#5252528F","#52525291","#52525294","#52525296","#52525299","#5252529C","#5252529E","#525252A1","#525252A3","#525252A6","#525252A8","#525252AD","#525252B0","#525252B3","#525252B5","#525252B8","#525252BA","#525252BD","#525252BF","#525252C2","#525252C4","#525252C7","#525252C9","#525252CC","#525252CF","#525252D1","#525252D4","#525252D6","#525252D9","#525252DB","#525252DE","#525252E0","#525252E3","#525252E6","#525252E8","#525252EB","#525252ED","#525252F0","#525252F2","#525252F5","#525252F7","#525252FA","#525252FC"), xlab = "", ylab = "", axes = FALSE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}

imagine4 <- function (comm, col = col, fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE, add = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = c("#00000000","#21908C03","#21908C05","#21908C08","#21908C0A","#21908C0D","#21908C0F","#21908C12","#21908C14","#21908C17","#21908C1A","#21908C1C","#21908C1F","#21908C21","#21908C24","#21908C26","#21908C29","#21908C2B","#21908C2E","#21908C30","#21908C33","#21908C36","#21908C38","#21908C3B","#21908C3D","#21908C40","#21908C42","#21908C45","#21908C47","#21908C4A","#21908C4D","#21908C4F","#21908C52","#21908C54","#21908C57","#21908C59","#21908C5C","#21908C5E","#21908C61","#21908C63","#21908C66","#21908C69","#21908C6B","#21908C6E","#21908C70","#21908C73","#21908C75","#21908C78","#21908C7A","#21908C7D","#21908C80","#21908C82","#21908C85","#21908C87","#21908C8A","#21908C8C","#21908C8F","#21908C91","#21908C94","#21908C96","#21908C99","#21908C9C","#21908C9E","#21908CA1","#21908CA3","#21908CA6","#21908CA8","#21908CAD","#21908CB0","#21908CB3","#21908CB5","#21908CB8","#21908CBA","#21908CBD","#21908CBF","#21908CC2","#21908CC4","#21908CC7","#21908CC9","#21908CCC","#21908CCF","#21908CD1","#21908CD4","#21908CD6","#21908CD9","#21908CDB","#21908CDE","#21908CE0","#21908CE3","#21908CE6","#21908CE8","#21908CEB","#21908CED","#21908CF0","#21908CF2","#21908CF5","#21908CF7","#21908CFA","#21908CFC"), xlab = "", ylab = "", axes = FALSE, add = TRUE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}

imagine5 <- function (comm, col = col, fill = TRUE, 
                      xlab = "", ylab = "", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE, add = TRUE) 
{
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(1, 6, 15, 1))
  
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), purp = c("#00000000","#44015403","#44015405","#44015408","#4401540A","#4401540D","#4401540F","#44015412","#44015414","#44015417","#4401541A","#4401541C","#4401541F","#44015421","#44015424","#44015426","#44015429","#4401542B","#4401542E","#44015430","#44015433","#44015436","#44015438","#4401543B","#4401543D","#44015440","#44015442","#44015445","#44015447","#4401544A","#4401544D","#4401544F","#44015452","#44015454","#44015457","#44015459","#4401545C","#4401545E","#44015461","#44015463","#44015466","#44015469","#4401546B","#4401546E","#44015470","#44015473","#44015475","#44015478","#4401547A","#4401547D","#44015480","#44015487","#4401548A","#4401548C","#4401548F","#44015491","#44015494","#44015496","#44015499","#4401549C","#4401549E","#440154A1","#440154A3","#440154A6","#440154A8","#440154AD","#440154B0","#440154B3","#440154B5","#440154B8","#440154BA","#440154BD","#440154BF","#440154C2","#440154C4","#440154C7","#440154C9","#440154CC","#440154CF","#440154D1","#440154D4","#440154D6","#440154D9","#440154DB","#440154DE","#440154E0","#440154E3","#440154E6","#440154E8","#440154EB","#440154ED","#440154F0","#440154F2","#440154F5","#440154F7","#440154FA","#440154FC"), xlab = "", ylab = "", axes = FALSE, add = TRUE)
  box()
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1, 
         cex.axis = 1, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2, 
         cex.axis = 1, lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}



imagine3(grays_only3)
imagine4(subblues)
imagine5(subpurples)

imagine3(log_grays_only1)
imagine4(log_blues_only1)
imagine5(log_purples_only1)

imagine3(subgreys)
imagine4(subblues)
imagine5(subpurples)

##### EXTRA FIGURES

#### PCA OF ENVIRONMENTAL VARIABLES (Not incl in paper) ---------------------------------
```{r, echo=FALSE}
### READ IN AND PREPARE DATA
wholeenv <- read.csv("wholeenv.csv", header = TRUE, stringsAsFactors = FALSE)
wholeenv <- wholeenv[wholeenv == 0] <- NA
row.names(wholeenv) <- wholeenv[,1]
wholeenv <- wholeenv[,-c(1)]
wholeenv$fetch <- log(wholeenv$fetch)
colnames(wholeenv) <- c("Chl a", "Dissolved O2", "Nitrates", "pH", "Phosphates", "Salinity", "Silicates", "SST", "Fetch", "Seagrass LAI", "Seagrass AFDM", "Algal AFDM", "Detritus AFDM") # just for vector arrows
pc <- princomp(wholeenv)
summary(pc) ## See proportion of variance explained by all principal axes... PC1 is 75.2% and PC2 is 14.4 %
pcasco <- scores(pc) #RECORD 1ST TWO PC SCORES
pca1 <- pcasco[,1] 
pca2 <- pcasco[,2] 
colrs <- vector() 

# SET SITE COLOURS ACCORDING TO PURPLE-BLUE-GREY COLOUR GROUPINGS FROM SPECIES CO-OCCURRENCE PLOT, MAKE PLOT

colrs <- c()
colrs[grep(c("HL"), row.names(pcasco))] <- "#21908C52" 
colrs[grep(c("IN"), row.names(pcasco))] <- "#21908C52" 
colrs[grep(c("RB"), row.names(pcasco))] <- "#21908C52" 
colrs[grep(c("RA"), row.names(pcasco))] <- "#525252F7"
colrs[grep(c("LH"), row.names(pcasco))] <- "#525252F7"
colrs[grep(c("SS"), row.names(pcasco))] <- "#525252F7"
colrs[grep(c("DK"), row.names(pcasco))] <- "#440154"
colrs[grep(c("EB"), row.names(pcasco))] <- "#440154"
colrs[grep(c("JB"), row.names(pcasco))] <- "#440154"
colrs[grep(c("GB"), row.names(pcasco))] <- "#440154"
colrs[grep(c("CB"), row.names(pcasco))] <- "#440154"
colrs[grep(c("DC"), row.names(pcasco))] <- "#440154"
colrs[grep(c("SA"), row.names(pcasco))] <- "#440154"
par(bty = "l") #set graphical parameters
plot.new()
plot.window(xlim = range(pca1) + c(-1,1), ylim = range(pca2) + c(-1,2)) #plot window
points(y=pca2,x=pca1,bg=colrs,pch =21, cex = 2) #plot points
text(y=pca2, x=pca1, labels = rownames(pcasco), cex = 0.7, pos = 3, col = "grey50") #label the site RDA points with the site names. 
axis(1, lwd = 0, lwd.ticks = 1) # set the x axis ticks and labels
axis(2, lwd =0, lwd.ticks =1) #set the y axis ticks and labels
title(xlab = "PCA1 (75.2%)", ylab = "PCA2 (14.4%)") #give axis titles
box() #close the graph box 
envfit <- envfit(pc, wholeenv) ##Add environmental vectors
plot(envfit, choices = c(1,2), arrow.mul = 6, at = c(0,0), axis = FALSE, 
     p.max = NULL, col = "black", cex = 0.8, add = TRUE)




#### PCA OF ENVIRONMENTAL VARIABLES (FIG 4B) ---------------------------------

### 1. READ IN AND PREPARE DATA
wholeenv <- read.csv("wholeenv.csv", header = TRUE, stringsAsFactors = FALSE)
wholeenv <- wholeenv[wholeenv == 0] <- NA
row.names(wholeenv) <- wholeenv[,1]
wholeenv <- wholeenv[,-c(1)]
wholeenv$fetch <- log(wholeenv$fetch)
colnames(wholeenv) <- c("Chl a", "Dissolved O2", "Nitrates", "pH", "Phosphates", "Salinity", "Silicates", "SST", "Fetch", "Seagrass LAI", "Seagrass AFDM", "Algal AFDM", "Detritus AFDM") # just for vector arrows

### 2. MAKE PCA OBJECT AND RECORD 1ST TWO PC SCORES

pc <- princomp(wholeenv)

summary(pc) ## See proportion of variance explained by all principal axes... PC1 is 75.2% and PC2 is 14.4 %! 
pcasco <- scores(pc) 
pca1 <- pcasco[,1] 
pca2 <- pcasco[,2] 
colrs <- vector() 


### 3. SET SITE COLOURS ACCORDING TO PURPLE-BLUE-GREY COLOUR GROUPINGS FROM SPECIES CO-OCCURRENCE PLOT

colrs <- c()
colrs[grep(c("HL"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("IN"), row.names(pcasco))] <- "#525252F7"
colrs[grep(c("RB"), row.names(pcasco))] <- "#440154" 
colrs[grep(c("RA"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("LH"), row.names(pcasco))] <- "#440154"
colrs[grep(c("SS"), row.names(pcasco))] <- "#440154"
colrs[grep(c("DK"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("EB"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("JB"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("GB"), row.names(pcasco))] <- "#440154"
colrs[grep(c("CB"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("DC"), row.names(pcasco))] <- "#21908C" 
colrs[grep(c("SA"), row.names(pcasco))] <- "#440154"

pch <- c()
pch[grep(c("HL"), row.names(pcasco))] <- 21
pch[grep(c("IN"), row.names(pcasco))] <- 23
pch[grep(c("RB"), row.names(pcasco))] <- 22
pch[grep(c("RA"), row.names(pcasco))] <- 21
pch[grep(c("LH"), row.names(pcasco))] <- 24
pch[grep(c("SS"), row.names(pcasco))] <- 24
pch[grep(c("DK"), row.names(pcasco))] <- 23
pch[grep(c("EB"), row.names(pcasco))] <- 23
pch[grep(c("JB"), row.names(pcasco))] <- 24
pch[grep(c("GB"), row.names(pcasco))] <- 24
pch[grep(c("CB"), row.names(pcasco))] <- 24
pch[grep(c("DC"), row.names(pcasco))] <- 22
pch[grep(c("SA"), row.names(pcasco))] <- 22

### 4. SET GRAPHICS PARAMETERS, PLOT POINTS ETC
dev.off()
par(bty = "l") 
plot.new()
plot.window(xlim = range(pca1) + c(-1,1), ylim = range(pca2) + c(-1,2))
points(y=pca2,x=pca1,bg=colrs,pch =pch, cex = 2) #plot the site rda points
text(y=pca2, x=pca1, labels = rownames(pcasco), cex = 0.7, pos = 3, col = "grey50") #label the site RDA points with the site names. 
axis(1, lwd = 0, lwd.ticks = 1) # set the x axis ticks and labels
axis(2, lwd =0, lwd.ticks =1) #set the y axis ticks and labels
title(xlab = "PCA1 (75.2%)", ylab = "PCA2 (14.4%)") #give axis titles
box() #close the graph box 


### 5. ADD ENV VECTORS
envfit <- envfit(pc, wholeenv)
plot(envfit, choices = c(1,2), arrow.mul = 6, at = c(0,0), axis = FALSE, 
     p.max = NULL, col = "black", cex = 0.8, add = TRUE)


# SUPPLEMENTAL: ggpairs -----------------------------------------------------------------

install.packages("GGally")
library(GGally)
X1 <- data.frame(X)
X1 <- X1[,c(4:7,15:24)]
ggpairs(X1)
par(mar = c(0,0,0,0))
plot.new()
color.legend(0.05,0.4,0.95,0.5, legend = c(0,0.25,0.5,0.75,1),cex = 2.5 ,align = "rb", rect.col = purp, gradient = "x")
color.legend(0.05,0.5,0.95,0.6,  rect.col = turq, gradient = "x")
color.legend(0.05,0.6,0.95,0.7,  rect.col = grs, gradient = "x")


preds <- predict(model)


## try new hmsc package
library(devtools)
install.packages('BayesLogit')
install_github('hmsc-r/HMSC')  

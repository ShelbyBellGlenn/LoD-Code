########################################################################################################
# Title: Estimating the uncertainty in deconvolution estimates for the extended library
# Topic: Analyses to assess uncertainty in deconvolution
# Authors: D. Koestler and S. Bell-Glenn
# Date: October 4th, 2022
########################################################################################################

########################################################################################################
# install and load necessary R packages
########################################################################################################
library(quadprog)
library(FlowSorted.BloodExtended.EPIC)
library(FlowSorted.Blood.EPIC)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(coxed)
library(bcaboot)
library(boot)

########################################################################################################
# read in the necessary source code and objects for deconvolution
########################################################################################################
source("IDOL revised functions 10_04_2018.R") # contains functions for deconvolution

########################################################################################################
# read in and prepare the necessary data files
########################################################################################################

FlowSorted.BloodExtended.EPIC
head(FlowSorted.BloodExtended.EPIC.compTable)

RGsetTargets <- FlowSorted.BloodExtended.EPIC[,FlowSorted.BloodExtended.EPIC$CellType == "MIX"]

sampleNames(RGsetTargets) <- paste(RGsetTargets$CellType,seq_len(dim(RGsetTargets)[2]), sep = "_")
RGsetTargets

#beta matrix corresponding to the cellular mixtures we have
betasmixtures<-getBeta(RGsetTargets)
dim(betasmixtures)#[1] 865859     12

#the true cell proportions
cellmix<-as.matrix(pData(RGsetTargets)[,colnames(FlowSorted.BloodExtended.EPIC.compTable)])

#get reference data to build library from
RGsetTargets2 <- FlowSorted.BloodExtended.EPIC[,FlowSorted.BloodExtended.EPIC$CellType != "MIX"]

betasreferences<-getBeta(RGsetTargets2)
colnames(betasreferences) <- RGsetTargets2$CellType
dim(betasreferences) #[1] 865859     56

#build NxP reference matrix where one col gives CellType
referencePd <- cbind(RGsetTargets2$Sample_Name, RGsetTargets2$Slide, 
                     RGsetTargets2$Array, RGsetTargets2$Basename, RGsetTargets2$Sample_ID, 
                     RGsetTargets2$CellType_long, RGsetTargets2$CellType,RGsetTargets2$sex)
rownames(referencePd) <- RGsetTargets2$Sample_Name
colnames(referencePd) <- c("Sample_Name", "Slide", "Array",
                           "Basename", "Sample_ID", "CellType_long",
                           "CellType", "sex")
referencePd <- data.frame(referencePd)

#the cells we have information for
cell_ident <- RGsetTargets2$CellType

#list of the 12 unique cell types
cells <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                     "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")


###############################################################################################
# Subset the data to the IDOL optimized CpGs only
###############################################################################################

load("IDOLOptimizedCpGsBloodExtended.rda")

#subset the cell specific reference data
extendend.idol.betas <- betasreferences[IDOLOptimizedCpGsBloodExtended,]

#subset the mixture reference data 
extended.idol.mix.betas <- betasmixtures[IDOLOptimizedCpGsBloodExtended,]


###############################################################################################
# Calculate bootstrap confidence intervals for each sample and cell type 
# METHOD: Basic Bootstrap Method
###############################################################################################

#number of bootstrap estimates to get
Nboot <- 10000

#matrix to save deconvolution the estimates 
sample.ests <- vector(mode="list", length=12)

for(s in 1:12){
  
  ests <- matrix(NA, nrow = Nboot, ncol = length(cells))
  colnames(ests) <- cells
  
  for(j in 1:Nboot){
  
    #get a random sample from each of the cell types to use as our matrix for deconvolution 
    random.cell.spec <- matrix(NA, nrow = length(IDOLOptimizedCpGsBloodExtended), ncol = length(cells))
    rownames(random.cell.spec) <- IDOLOptimizedCpGsBloodExtended
    colnames(random.cell.spec) <- cells
    for(i in 1:length(cells)){
      
      #get the cols from the reference matrix corresponding to cell type i
      tmp <- extendend.idol.betas[,colnames(extendend.idol.betas) == cells[i]]
      
      #get the number of cols in this matrix and then select a random sample to use
      k <- ncol(tmp)
      rand.col <- sample(1:k, 1)
      
      #save this sample in our matrix to use for deconvolution later
      random.cell.spec[,i] <- tmp[,rand.col]
    }
    
    #deconvolute one of the reconstructed mixtures using this matrix
    cellpredictions = projectWBCnew(as.matrix(extended.idol.mix.betas[,s]), random.cell.spec)
    
    
    #save the deconvolution estimates as percentages
    ests[j,] <- cellpredictions*100
  }
  sample.ests[[s]] <- ests
}


#calculate the bootsrap estimated variance for each reconstructed mixture and each cell type
boot.est.var <- matrix(NA, nrow = 12, ncol = 12)
rownames(boot.est.var) <- colnames(extended.idol.mix.betas)
colnames(boot.est.var) <- cells
for(i in 1:12){
  tmp <- sample.ests[[i]]/100
  for(j in 1:12){
    v <- var(tmp[,j])
    boot.est.var[i,j] <- v
  }
}

save(boot.est.var, file = "Bootstrap Est Var of Deconvolution Ests.RData")

#order the columns of the matrices for each cell type and get the upper and lower limit for CIs
sample.ests.ordered <- vector(mode="list", length=12)
for(i in 1:12){
  sample.ests.ordered[[i]] <- apply(sample.ests[[i]], 2, sort)
}


# Load the LoD estimates using concentration parameter 73 and cut the lower limit of 
# confidence intervals off at the LoD
#load("LoDResults_73.RData")

#NOTE: LoD ests are %

#matrix to store confidence intervals
# cis <- matrix(NA, nrow=12, ncol = 12)
# colnames(cis) <- cells
# rownames(cis) <- colnames(extended.idol.mix.betas)
# 
# for(k in 1:12){
#   c <- NULL
#   for(i in 1:length(cells)){
#     
#     #get the lower limit and cut off at LoD
#     lower <- sample.ests.ordered[[k]][,i][250]
#     if(lower < LoDResults_73[i]){
#       lower <- LoDResults_73[i]
#     }
#     
#     #get the upper limit
#     upper <- sample.ests.ordered[[k]][,i][9750]
#     
#     #save the estimate in a easily readable format
#     c[i] <- paste("(", round(lower,4), ", ", round(upper,4), ")", sep="")
#   }
#   cis[k,] <- c
# }

#save these cis
#save(cis, file = "Bootstrap CIs cut off.RData")

#now calculate but dont cut off at LoD

#matrix to store confidence intervals
cis2 <- matrix(NA, nrow=12, ncol = 12)
colnames(cis2) <- cells
rownames(cis2) <- colnames(extended.idol.mix.betas)

for(k in 1:12){
  c <- NULL
  for(i in 1:length(cells)){
    
    #get the lower 
    lower <- sample.ests.ordered[[k]][,i][250]
    
    #get the upper limit
    upper <- sample.ests.ordered[[k]][,i][9750]
    
    #save the estimate in a easily readable format
    c[i] <- paste("(", round(lower,4), ", ", round(upper,4), ")", sep="")
  }
  cis2[k,] <- c
}

#save these estimates too
#save(cis2, file = "Bootstrap CIs.RData")
load("Bootstrap CIs.RData")

# Calculate coverage of our bootstrap confidence intervals using true props

covered <- rep(1, 144)
indx <- 1
for(k in 1:12){
  
  for(i in 1:length(cells)){
    
    #get the lower limit and cut off at LoD
    lower <- sample.ests.ordered[[k]][,i][250]
    
    #get the upper limit
    upper <- sample.ests.ordered[[k]][,i][9750]
    
    #get the true value for this mixture, for this cell type
    true.value <- cellmix[k,i]*100
    
    #check if true value is contained within the lower and upper limits
    if(true.value < lower || true.value > upper ){
      covered[indx] <- 0
    }
    indx <- indx + 1 
  }
  
}

#now calculate coverage %
coverage <- sum(covered)/144

# Calculate coverage of our bootstrap confidence intervals using decon ests

load("ExtendedMeanMethylation.RData")
load("Deconvolution Estimates Extended.RData")

get_cis <- function(vec){
  
  lower <- NULL
  upper <- NULL
  for(i in 1:12){
    tmp <- strsplit(vec[i], split = ",", fixed = T)
    
    tmp1 <- tmp[[1]][1]
    lower[i] <- as.numeric(gsub("[(]", "", tmp1))
    
    tmp2 <- tmp[[1]][2]
    upper[i] <- as.numeric(gsub("[)]", "", tmp2))
    
  }
  
  return(list(lo = lower, up = upper))
  
}

#create object to calculate coverage

covered <- rep(1, 144)
indx <- 1
for(k in 1:12){
  
  c <- get_cis(cis2[k,])
  
  lower <- c$lo
  upper <- c$up  
    
  for(i in 1:12){
    est.value <- round(decon.statistic[k,i], 5)
    
    if(est.value < lower[i] || est.value > upper[i] ){
      covered[indx] <- 0
    }
    indx <- indx + 1 
  }
}

sum(covered)/144

##############################################################
# Calculate the average coverage of the confidence intervals
#
##############################################################


get_cis <- function(vec){
  
  lower <- NULL
  upper <- NULL
  for(i in 1:12){
    tmp <- strsplit(vec[i], split = ",", fixed = T)
    
    tmp1 <- tmp[[1]][1]
    lower[i] <- as.numeric(gsub("[(]", "", tmp1))
    
    tmp2 <- tmp[[1]][2]
    upper[i] <- as.numeric(gsub("[)]", "", tmp2))
    
  }
  
  return(list(lo = lower, up = upper))
  
}

width <- NULL
for(i in 1:12){
  
  c <- get_cis(cis2[i,])
  
  lower <- c$lo
  upper <- c$up
 
  diff <- upper-lower
  
  width[i] <- mean(diff)
}



###############################################################################################
# Calculate bootstrap confidence intervals for each sample and cell type 
# METHOD: Semi Parametric Bootstrap Method
###############################################################################################

#number of bootstrap estimates to get
Nboot <- 10000

# download the beta parameter estimates for the 1200 idol cpgs
load("UpdatedExtenedIDOLCellSpecParams.RData")

#initialize matrix to save cell predictions
sample.ests <- vector(mode="list", length=12)
for (i in 1:12) {
  sample.ests[[i]] <- matrix(NA, nrow = Nboot, ncol = length(cells))
  colnames(sample.ests[[i]]) <- cells
}


for(j in 1:Nboot){
  
  #initialize ransom matrix to use for deconvolution
  random.cell.spec <- matrix(NA, nrow = length(IDOLOptimizedCpGsBloodExtended), ncol = length(cells))
  rownames(random.cell.spec) <- IDOLOptimizedCpGsBloodExtended
  colnames(random.cell.spec) <- cells
  
  #generate a random matrix to use for deconvolution 
  for(i in 1:length(cells)){
      tmp = apply(idol.params[[i]], 1, function(p) rbeta(1, shape1 = p[1], shape2 = p[2]))
      random.cell.spec[,i] = tmp
  }
  
  cellpredictions = projectWBCnew(extended.idol.mix.betas, random.cell.spec)*100
  
  #save the deconvolution estimates as percentages in our list for each reconstructed mixture
  for(k in 1:12){
    sample.ests[[k]][j,] <- cellpredictions[k,]
  }

}

# Load the LoD estimates using concentration parameter 73 and cut the lower limit of 
# confidence intervals off at the LoD
load("LoDResults_73.RData")

#order the columns of the matrices for each cell type and get the upper and lower limit for CIs
sample.ests.ordered <- vector(mode="list", length=12)
for(i in 1:12){
  sample.ests.ordered[[i]] <- apply(sample.ests[[i]], 2, sort)
}

#NOTE: LoD ests are %

#matrix to store confidence intervals
cis <- matrix(NA, nrow=12, ncol = 12)
colnames(cis) <- cells
rownames(cis) <- colnames(extended.idol.mix.betas)

for(k in 1:12){
  c <- NULL
  for(i in 1:length(cells)){
    
    #get the lower limit and cut off at LoD
    lower <- sample.ests.ordered[[k]][,i][250]
    if(lower < LoDResults_73[i]){
      lower <- LoDResults_73[i]
    }
    
    #get the upper limit
    upper <- sample.ests.ordered[[k]][,i][9750]
    
    #save the estimate in a easily readable format
    c[i] <- paste("(", round(lower,4), ", ", round(upper,4), ")", sep="")
  }
  cis[k,] <- c
}

#save these cis
#save(cis, file = "Bootstrap Semi-Parametric CIs cut off.RData")

#now calculate but dont cut off at LoD

#matrix to store confidence intervals
cis2 <- matrix(NA, nrow=12, ncol = 12)
colnames(cis2) <- cells
rownames(cis2) <- colnames(extended.idol.mix.betas)

for(k in 1:12){
  c <- NULL
  for(i in 1:length(cells)){
    
    #get the lower limit and cut off at LoD
    lower <- sample.ests.ordered[[k]][,i][250]
    
    #get the upper limit
    upper <- sample.ests.ordered[[k]][,i][9750]
    
    #save the estimate in a easily readable format
    c[i] <- paste("(", round(lower,4), ", ", round(upper,4), ")", sep="")
  }
  cis2[k,] <- c
}

#save these estimates too
#save(cis2, file = "Bootstrap Semi-Parametric CIs.RData")

# Calculate coverage of our bootstrap confidence intervals

covered <- rep(1, 144)
indx <- 1
for(k in 1:12){
  
  for(i in 1:length(cells)){
    
    #get the lower limit and cut off at LoD
    lower <- sample.ests.ordered[[k]][,i][250]
    
    #get the upper limit
    upper <- sample.ests.ordered[[k]][,i][9750]
    
    #get the true value for this mixture, for this cell type
    true.value <- cellmix[k,i]*100
    
    #check if true value is contained within the lower and upper limits
    if(true.value < lower || true.value > upper ){
      covered[indx] <- 0
    }
    indx <- indx + 1 
  }
  
}

#now calculate coverage %
coverage <- sum(covered)/144







###############################################################################################
# Check whether or not our deconvolution estimates are below the LoD/LoB for the mixtures 
# that are truly missing a cell type
###############################################################################################


#load mean methylation for the extended library
load("ExtendedMeanMethylation.RData")

#load limit of blank ests
load("LoB Results New Con 73.RData")


#deconvolute the reconstructed mixtures
cellpredictionscheck = projectWBCnew(extended.idol.mix.betas, meanmeth)*100

rounded.ests <- round(cellpredictionscheck, 5)

#look at the predictions, the true values and the LoDs
rounded.ests
cellmix*100
LoBResults_73
LoDResults_73

###############################################################################################
# Calculate bootstrap confidence intervals for each sample and cell type 
# METHOD: Bias Corrected and accelerated Bootstrap Method
###############################################################################################


###########
# Step 1 
# Get the original "statistic" which is the cell proportion estimates using the mean methylation of the idol cpgs
###########

#load("ExtendedMeanMethylation.RData")
#decon.statistic <- projectWBCnew(extended.idol.mix.betas, meanmeth)*100
#save(decon.statistic, file = "Deconvolution Estimates Extended.RData")


###########
# Step 2
# Get our bootstrapped estimates of the statistic
#
# This was done above using the normal non parametric bootstrap procedure
###########



###########
# Step 3
# Estimate the acceleration parameter, which is proportional to the skewness of the bootstrap distribution
###########



###########
# Step 4
# Estimate the bias-corrected parameter, which is related to the proportion of bootstrap estimates less than 
# the observed statistic
###########

# We do this for each sample, each cell type. This will give us a 12x12 matrix of bias corrected 




###########
# Step 5
# Get the adjusted CIs based on the adjusted quantiles
###########






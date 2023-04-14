######################################################################################
# Recalculate necesssary parameters to calculate LoB and LoD to 
######################################################################################


#####################################################################
# Read in necessary packages
####################################################################

library(quadprog)
library(betareg)
library(simstudy)
library(rBeta2009)
library(FlowSorted.BloodExtended.EPIC)
library(FlowSorted.Blood.EPIC)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(EnvStats)
library(ggplot2)
library(doParallel)

###################################################################
#read in and prepare full extended data
###################################################################
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
uniquecelltypes <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                     "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")


###################################################################
# Load in mean and variance matrices
###################################################################

load('FullExtendedMeanMeth.RData')
load('FullExtendedVarMeth.RData')

###################################################################
# Subset reference data and mixture data to only the idol cpgs
###################################################################

#load the 1200 idol cpgs
load("IDOLOptimizedCpGsBloodExtended.rda")

#subset the cell specific reference data
extendend.idol.betas <- betasreferences[IDOLOptimizedCpGsBloodExtended,]

#subset the mixture reference data 
extended.idol.mix.betas <- betasmixtures[IDOLOptimizedCpGsBloodExtended,]

#subset mean and variance methylation to only the idol cpgs
idol.mean.meth <- meanmeth[IDOLOptimizedCpGsBloodExtended,]
idol.var.meth <- varmeth[IDOLOptimizedCpGsBloodExtended,]


###################################################################
# Estimate shape parameters for beta distribution via MME
###################################################################

#get parameter estimates for the idol cpgs
idol.params <- vector(mode="list", length=length(uniquecelltypes))
for(i in 1:length(idol.params)){
  alphas <- NULL
  betas <- NULL
  for(j in 1:nrow(extendend.idol.betas)){
    params <- try(ebeta(extendend.idol.betas[j, colnames(extendend.idol.betas) == uniquecelltypes[i]], 
                        method = "mme"), silent = T )
    if( !(is(params) == "try-error")  ){
      alphas[j] <- params$parameters[1]
      betas[j] <- params$parameters[2]
    }
    else{
      alphas[j] = betas[j] = NA
    }
  }
  
  temp <- cbind(alphas, betas)
  rownames(temp) <- rownames(extendend.idol.betas)
  colnames(temp) <- c("alphas","betas")
  
  #return the matrix
  idol.params[[i]] <- temp
}


#save these estimates
names(idol.params) <- c("Bas", "Bmem", "Bnv","CD4mem", 
                        "CD4nv" ,"CD8mem", "CD8nv", 
                        "Eos", "Mono", "Neu", 
                        "NK", "Treg")


###################################################################
# Check to see which parameters have bad estimates for each cell type
###################################################################

bad.shape.ests <- vector(mode="list", length=length(uniquecelltypes))
for(i in 1:length(idol.params)){
  ind2 <- 1
  for(j in 1:nrow(extendend.idol.betas)){
    if(idol.params[[i]][j,1] >= 100000 || idol.params[[i]][j,2] >= 100000 ){
      bad.shape.ests[[i]][ind2] <- rownames(idol.params[[i]])[j]
      ind2 <- ind2 + 1
    }
  }
}


###################################################################
# For the idol cpgs with bad shape parameter estimates in a cell type, we will 
# add a small amount of normal noise to the cell types that the issues
# arise in, and then recalculate the shape parameters
###################################################################


for (i in 1:length(bad.shape.ests)) {
  indx <- which(rownames(extendend.idol.betas) %in% bad.shape.ests[[i]])
  if(length(bad.shape.ests[[i]]) > 0 ){
    for(j in 1:length(bad.shape.ests[[i]])){
      cols <- colnames(extendend.idol.betas) == uniquecelltypes[i]
      for(s in 1:length(extendend.idol.betas[indx[j],cols])){
        extendend.idol.betas[indx[j],cols][s] <- extendend.idol.betas[indx[j],cols][s] + rnorm(n=1, mean = 0, sd=.01)
        #check to make sure we don't have a values outside of the 0-1 range
        if(extendend.idol.betas[indx[j],cols][s] >= 1){
          extendend.idol.betas[indx[j],cols][s] <- 0.99
        }
        if(extendend.idol.betas[indx[j],cols][s] <= 0){
          extendend.idol.betas[indx[j],cols][s] <- 0.001
        }
      }
    }
  }
}

#get parameter estimates for the idol cpgs
idol.params <- vector(mode="list", length=length(uniquecelltypes))
for(i in 1:length(idol.params)){
  alphas <- NULL
  betas <- NULL
  for(j in 1:nrow(extendend.idol.betas)){
    params <- try(ebeta(extendend.idol.betas[j, colnames(extendend.idol.betas) == uniquecelltypes[i]], 
                        method = "mme"), silent = T )
    if( !(is(params) == "try-error")  ){
      alphas[j] <- params$parameters[1]
      betas[j] <- params$parameters[2]
    }
    else{
      alphas[j] = betas[j] = NA
    }
  }
  
  temp <- cbind(alphas, betas)
  rownames(temp) <- rownames(extendend.idol.betas)
  colnames(temp) <- c("alphas","betas")
  
  #return the matrix
  idol.params[[i]] <- temp
}


#save these estimates
names(idol.params) <- c("Bas", "Bmem", "Bnv","CD4mem", 
                        "CD4nv" ,"CD8mem", "CD8nv", 
                        "Eos", "Mono", "Neu", 
                        "NK", "Treg")

#save these updated parameters
#save(idol.params, file = "UpdatedExtenedIDOLCellSpecParams.RData")


###################################################################
# Try to simulate cell specific methylation data using these 
# estimates and check to see if we have the behavior we expect 
###################################################################


N = 10 # number of in-silico mixtures
cells = uniquecelltypes # cell type names
K = length(cells) # number of cell types
J = length(IDOLOptimizedCpGsBloodExtended) # number of CpGs 
concentration.param = c(18,73,127) # Dirichlet concentration parameter as per Meier et al., (2021)
true.w = matrix(NA, nrow = N, ncol = K) # matrix to store "true" cell proportions for the in-silico mixtures
colnames(true.w) = cells
rownames(true.w) = paste("sample_", 1:N, sep = "")

# 2. generate cell-specific methylation data for each in-silico sample
cellspecific.methylation = vector(mode="list", length=N)
for(i in 1:N) {
  tmp = lapply(idol.params, function(w)
    apply(w, 1, function(p) rbeta(1, shape1 = p[1], shape2 = p[2])))
  cellspecific.methylation[[i]] = do.call(cbind, tmp)
}
names(cellspecific.methylation) = rownames(true.w)	   



#the simulated data is looking better for the samples I think!!


###################################################################
# Estimate precision parameter estimates using mixture data for the
# idol cpgs
###################################################################


cellmix <- as.data.frame(cellmix)

Phi.ests <- NULL
for(i in 1:nrow(extended.idol.mix.betas)){
 mod <- betareg(extended.idol.mix.betas[i,]~ cellmix$Bas + cellmix$Bmem + cellmix$Bnv +
                 cellmix$CD4mem + cellmix$CD4nv + cellmix$CD8mem + cellmix$CD4nv + 
                 cellmix$Eos + cellmix$Mono + cellmix$Neu + cellmix$NK,
               data= as.data.frame(extended.idol.mix.betas))
Phi.ests[i] <- mod$coefficients$precision
}

names(Phi.ests) <- rownames(extended.idol.mix.betas)


###################################################################
# Check to see which cpgs have bad phi estimates
###################################################################


bad.phi.ests <- NULL
ind <- 1
for(i in 1:length(Phi.ests)){
  if(Phi.ests[i] >= 100000 || is.na(Phi.ests[i]) ){
    bad.phi.ests[ind] <- names(Phi.ests)[i]
    ind <- ind + 1
  }
}



###################################################################
# Replace the 282 bad precision estimates by sampling with replacement
# from the 918 valid phi estimates 
###################################################################

#which phi estimates are okay?
valid.phis <- Phi.ests[which(Phi.ests < 100000)]
valid.phis <- unname(valid.phis)

#get the index of the cpgs that have bad precision parameter estimates so we can replace them
indx2 <- which(rownames(extendend.idol.betas) %in% bad.phi.ests)

#replace the bad estimates with a random good one
for(i in 1:length(bad.phi.ests)){
  Phi.ests[indx2[i]] <- sample(valid.phis, 1, replace = T)
}

summary(Phi.ests)

#save the updated phi estimates
#save(Phi.ests, file = "UpdatedValidPhiEsts.RData")

###################################################################
# Try to simulate mixture methylation data using these estimates 
# and check to see if we have the behavior we expect 
###################################################################

#generate "true" cell proportions from a dirichlet distribution
dirichlet.params = runif(K) 
dirichlet.params = dirichlet.params/sum(dirichlet.params)*concentration.param[2]
true.w[,] = rdirichlet(N, shape = dirichlet.params)

# generate methylation data for in-silico mixtures
mixture.data = matrix(NA, nrow = J, ncol = N)
colnames(mixture.data) = rownames(true.w)
rownames(mixture.data) = IDOLOptimizedCpGsBloodExtended
for(i in 1:N) {
  mean.meth = cellspecific.methylation[[i]] %*% true.w[i,]
  beta.shape.parameters = do.call(cbind, betaGetShapes(as.vector(mean.meth), Phi.ests))
  mixture.data[,i] = apply(beta.shape.parameters, 1, function(w) rbeta(1, shape1 = w[1], shape2 = w[2])) 
}



########################################################################################################
# Title: Estimating the limit of blank for DNA methylation-based cell mixture deconvolution
#        for the extended library
# Topic: Analyses to assess the detection limit
# Authors: D. Koestler and S. Bell-Glenn
# Date: March 17th, 2022
########################################################################################################

########################################################################################################
# install and load necessary R packages
########################################################################################################
library(quadprog)
library(betareg)
library(simstudy)
library(rBeta2009)
library(FlowSorted.BloodExtended.EPIC)
library(FlowSorted.Blood.EPIC)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)

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
uniquecelltypes <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv",
                     "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")

###############################################################################################
# Subset the data to the IDOL optimized CpGs only
##############################################################################################

load("IDOLOptimizedCpGsBloodExtended.rda")

#subset the cell specific reference data
extendend.idol.betas <- betasreferences[IDOLOptimizedCpGsBloodExtended,]

#subset the mixture reference data 
extended.idol.mix.betas <- betasmixtures[IDOLOptimizedCpGsBloodExtended,]

#load mean methylation information for the 1200 idol cpgs
load("ExtendedMeanMethylation.RData")

########################################################################################################
# estimate the shape parameters for the beta distribution using mean and variance
########################################################################################################

#load("ExtenedIDOLCellSpecParams.RData")


load("UpdatedExtenedIDOLCellSpecParams.RData")

########################################################################################################
# estimate the precisions for each of the IDOL CpGs for the extended library
########################################################################################################
# cellmix <- as.data.frame(cellmix)
# 
# Phi.ests <- NULL
# for(i in 1:nrow(extended.idol.mix.betas)){
#   mod <- betareg(extended.idol.mix.betas[i,]~ cellmix$Bas + cellmix$Bmem + cellmix$Bnv +
#                   cellmix$CD4mem + cellmix$CD4nv + cellmix$CD8mem + cellmix$CD4nv + 
#                   cellmix$Eos + cellmix$Mono + cellmix$Neu + cellmix$NK,
#                 data= as.data.frame(extended.idol.mix.betas))
#  Phi.ests[i] <- mod$coefficients$precision
# }

# names(Phi.ests) <- rownames(extended.idol.mix.betas)
#save so we don't have to run this code again
#save(Phi.ests, file = "ExtendedIDOLPhiests.Rdata")


# load("ExtendedIDOLPhiests.Rdata")
# #since 200 of the phi ests for the idol cpgs are too big... we will sample 
# #with replacement to get 1200 valid values to use to calculate the LoB
# valid.phis <- Phi.ests[which(Phi.ests < 100000)]
# valid.idol.phis <-NULL
# for(i in 1:nrow(extended.idol.mix.betas)){
#   valid.idol.phis[i] <- sample(valid.phis, 1, replace = T)
# }
# 
# #save these random valid phi ests to be used for the limit of detection
# save(valid.idol.phis, file = "RandomValidPhiEsts.RData")

#load("RandomValidPhiEsts.RData")

load("UpdatedValidPhiEsts.RData")

########################################################################################################
# simulation study to estimate the detection limit of cell mixture deconvolution
########################################################################################################

# define pertinent parameters for simulation study
N = 100 # number of in-silico mixtures
cells = uniquecelltypes # cell type names
K = length(cells) # number of cell types
J = length(IDOLOptimizedCpGsBloodExtended) # number of CpGs 
concentration.param = c(18,73,127) # Dirichlet concentration parameter as per Meier et al., (2021)
true.w = matrix(NA, nrow = N, ncol = K) # matrix to store "true" cell proportions for the in-silico mixtures
colnames(true.w) = cells
rownames(true.w) = paste("sample_", 1:N, sep = "")

LoBResults_127 = matrix(NA, nrow = 1, ncol = length(cells))
rownames(LoBResults_127) = "Conc_18"
colnames(LoBResults_127) = cells

for(k in 1:K) {
  true.w[,cells[k]] = 0  # set the "true" proportion of the ghost cell type (GCT) equal to zero
  
  # 1. generate "true" cell proportions from a dirichlet distribution for the non-GCTs
  dirichlet.params = runif(K-1) 
  dirichlet.params = dirichlet.params/sum(dirichlet.params)*concentration.param[3]
  true.w[,!(cells == cells[k])] = rdirichlet(N, shape = dirichlet.params)
  print(paste("Done with step 1, cell = ", k, sep = ""))
  
  # 2. generate cell-specific methylation data for each in-silico sample
  cellspecific.methylation = vector(mode="list", length=N)
  for(i in 1:N) {
    tmp = lapply(idol.params, function(w)
      apply(w, 1, function(p) rbeta(1, shape1 = p[1], shape2 = p[2])))
    cellspecific.methylation[[i]] = do.call(cbind, tmp)
  }
  names(cellspecific.methylation) = rownames(true.w)	    
  print(paste("Done with step 2, cell = ", k, sep = ""))
  
  # 3. generate methylation data for in-silico mixtures
  mixture.data = matrix(NA, nrow = J, ncol = N)
  colnames(mixture.data) = rownames(true.w)
  rownames(mixture.data) = IDOLOptimizedCpGsBloodExtended
  for(i in 1:N) {
    mean.meth = cellspecific.methylation[[i]] %*% true.w[i,]
    beta.shape.parameters = do.call(cbind, betaGetShapes(as.vector(mean.meth), Phi.ests))
    mixture.data[,i] = apply(beta.shape.parameters, 1, function(w) rbeta(1, shape1 = w[1], shape2 = w[2])) 
  }
  print(paste("Done with step 3, cell = ", k, sep = ""))
  
  # 4. deconvolute the in-silico mixtures and record performance
  cellpredictions = projectWBCnew(mixture.data, meanmeth[,cells])
  print(paste("Done with step 4, cell = ", k, sep = ""))
  
  # save the true cell proportions and predicted cell proportions for later use
  filename = paste("GCT_", cells[k], ".RData", sep = "")
  save(cellpredictions, true.w, file = filename)
  
  #calculate and save the limit of blank for this cell time
  LoBResults_127[k] = as.numeric(sort(cellpredictions[,cells[k]]*100)[95])
}


save(LoBResults_127, file = "LoB Results New Con 127.RData")

########################################################################################################
# visualize/analyze the results for concentration parameter of 73
########################################################################################################


#list of files with gct information
files = c("GCT_Bas.RData","GCT_Bmem.RData", "GCT_Bnv.RData","GCT_CD4mem.RData",
          "GCT_CD4nv.RData","GCT_CD8mem.RData","GCT_CD8nv.RData","GCT_Eos.RData",
          "GCT_Mono.RData","GCT_Neu.RData","GCT_NK.RData","GCT_Treg.RData")
name.file = as.character(sapply(files, function(w) strsplit(w, "[.]")[[1]][1]))
cell.name = as.character(sapply(name.file, function(w) strsplit(w, "_")[[1]][2]))

par(mfrow = c(3,4))

#bas
load(files[1])
hist(cellpredictions[,cell.name[1]]*100, xlab = "Cell Percentage (%)", 
    main = name.file[1], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Bas"], lty = "dotted", col = "red")
text(1 , 6, labels = paste("LoB = ", round(LoBResults_127[,"Bas"],3), sep = ""), 
       cex = 1.5, col = "red")
#dev.off()

#bmem
load(files[2])
hist(cellpredictions[,cell.name[2]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[2], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Bmem"], lty = "dotted", col = "red")
text(1 , 6, labels = paste("LoB = ", round(LoBResults_127[,"Bmem"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#bnv
load(files[3])
hist(cellpredictions[,cell.name[3]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[3], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Bnv"], lty = "dotted", col = "red")
text(1, 6, labels = paste("LoB = ", round(LoBResults_127[,"Bnv"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#cd4mem
load(files[4])
hist(cellpredictions[,cell.name[4]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[4], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"CD4mem"], lty = "dotted", col = "red")
text(1.25, 4, labels = paste("LoB = ", round(LoBResults_127[,"CD4mem"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#cd4nv
load(files[5])
hist(cellpredictions[,cell.name[5]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[5], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"CD4nv"], lty = "dotted", col = "red")
text(1, 2, labels = paste("LoB = ", round(LoBResults_127[,"CD4nv"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#cd8mem
load(files[6])
hist(cellpredictions[,cell.name[6]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[6], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"CD8mem"], lty = "dotted", col = "red")
text(1.25, 4, labels = paste("LoB = ", round(LoBResults_127[,"CD8mem"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#cd8nv
load(files[7])
hist(cellpredictions[,cell.name[7]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[7], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"CD8nv"], lty = "dotted", col = "red")
text(1.25, 4, labels = paste("LoB = ", round(LoBResults_127[,"CD8nv"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#eos
load(files[8])
hist(cellpredictions[,cell.name[8]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[8], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Eos"], lty = "dotted", col = "red")
text(1, 6, labels = paste("LoB = ", round(LoBResults_127[,"Eos"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#mono
load(files[9])
hist(cellpredictions[,cell.name[9]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[9], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Mono"], lty = "dotted", col = "red")
text(1, 6, labels = paste("LoB = ", round(LoBResults_127[,"Mono"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#neu
load(files[10])
hist(cellpredictions[,cell.name[10]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[10], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Neu"], lty = "dotted", col = "red")
text(1, 6, labels = paste("LoB = ", round(LoBResults_127[,"Neu"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#nk
load(files[11])
hist(cellpredictions[,cell.name[11]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[11], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"NK"], lty = "dotted", col = "red")
text(1, 6, labels = paste("LoB = ", round(LoBResults_127[,"NK"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()

#treg
load(files[12])
hist(cellpredictions[,cell.name[12]]*100, xlab = "Cell Percentage (%)", 
     main = name.file[12], xlim = c(0,2), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
abline(v = LoBResults_127[,"Treg"], lty = "dotted", col = "red")
text(1.25, 4, labels = paste("LoB = ", round(LoBResults_127[,"Treg"],3), sep = ""), 
     cex = 1.5, col = "red")
#dev.off()



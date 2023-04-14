########################################################################################################
# Title: Estimating the detection limit of DNA methylation-based cell mixture deconvolution
#        for the extended library
# Topic: Analyses to assess the detection limit
# Authors: D. Koestler and S. Bell-Glenn
# Date: March 17th, 2022
########################################################################################################

#NOTE: The following was done for a concentration parameter of 73

########################################################################################################
# install and load necessary R packages
########################################################################################################
library(quadprog)
library(simstudy)
library(rBeta2009)
library(FlowSorted.BloodExtended.EPIC)
library(FlowSorted.Blood.EPIC)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)

########################################################################################################
# read in the necessary data files
########################################################################################################

#different LoB results
load("LoB Results New Con 18.RData") # read in the table of the limits of blanks
load("LoB Results New Con 73.RData") # read in the table of the limits of blanks
load("LoB Results New Con 127.RData") # read in the table of the limits of blanks

load("UpdatedValidPhiEsts.RData")# precision estimates used to calculate the LOB
load("UpdatedExtenedIDOLCellSpecParams.RData") #load valid cell specific beta shape params for the 1200 idol cpgs
load("ExtendedMeanMethylation.RData") #load mean methylation information for the 1200 idol cpgs
load("IDOLOptimizedCpGsBloodExtended.rda") #load 1200 idol cpg names

#read in and prepare full extended data
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

#now subset to only the IDOL cpgs

#subset the cell specific reference data
extendend.idol.betas <- betasreferences[IDOLOptimizedCpGsBloodExtended,]

#subset the mixture reference data 
extended.idol.mix.betas <- betasmixtures[IDOLOptimizedCpGsBloodExtended,]


########################################################################################################
# read in the necessary source code and objects for deconvolution
########################################################################################################
source("IDOL revised functions 10_04_2018.R") # contains functions for deconvolution

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

#the type 2 error we want to control for and how much we will increment
#at each iteration until we find the LoD
t2.error = 0.20
increment = 0.0005

#create matrix to store results
LoDResults_127 = matrix(NA, nrow = 1, ncol = length(cells))
rownames(LoDResults_127) = "Conc_73"
colnames(LoDResults_127) = cells

cellpredictions_GCT = list()
truecellprops_GCT = list()

for(k in 1:K) {
  
  prop.GST = LoBResults_127[1, cells[k]]/100 + increment
  beta.t2 = 1.0
  
  while(beta.t2 > t2.error) {
    true.w[,cells[k]] = prop.GST  # set the "true" proportion of the ghost cell type (GCT) to be a small value
    
    # 1. generate "true" cell proportions from a dirichlet distribution for the non-GCTs
    dirichlet.params = runif(K-1) 
    
    dirichlet.params = dirichlet.params/sum(dirichlet.params)*concentration.param[2]
    cell.props.non.GST = rdirichlet(N, shape = dirichlet.params)
    cell.props.non.GST.scaled = (1-true.w[,cells[k]])*cell.props.non.GST
    true.w[,!(cells == cells[k])] = cell.props.non.GST.scaled
    
    #print(paste("Done with step 1, cell = ", k, sep = ""))
    
    # 2. generate cell-specific methylation data for each in-silico sample
    cellspecific.methylation = list()
    for(i in 1:N) {
      tmp = lapply(idol.params, function(w)
        apply(w, 1, function(p) rbeta(1, shape1 = p[1], shape2 = p[2])))
      cellspecific.methylation[[i]] = do.call(cbind, tmp)
    }
    names(cellspecific.methylation) = rownames(true.w)	    
    #print(paste("Done with step 2, cell = ", k, sep = ""))
    
    # 3. generate methylation data for in-silico mixtures
    mixture.data = matrix(NA, nrow = J, ncol = N)
    colnames(mixture.data) = rownames(true.w)
    rownames(mixture.data) = rownames(extended.idol.mix.betas)
    for(i in 1:N) {
      mean.meth = cellspecific.methylation[[i]] %*% true.w[i,]
      beta.shape.parameters = do.call(cbind, betaGetShapes(as.vector(mean.meth), Phi.ests))
      mixture.data[,i] = apply(beta.shape.parameters, 1, function(w) rbeta(1, shape1 = w[1], shape2 = w[2])) 
    }
    #print(paste("Done with step 3, cell = ", k, sep = ""))
    
    # 4. deconvolute the in-silico mixtures and record performance
    cellpredictions = projectWBCnew(mixture.data, meanmeth[,cells])*100
    #print(paste("Done with step 4, cell = ", k, sep = ""))
    
    # 5. compute the beta parameter (type-2-error)
    LoB = LoBResults_127[1, cells[k]]
    beta.t2 = mean(cellpredictions[,cells[k]]<= LoB)
    prop.GST = prop.GST + increment	
  }
  LoD = (prop.GST-increment)*100
  LoDResults_127[1,cells[k]] = LoD
  cellpredictions_GCT[[k]] = cellpredictions/100
  truecellprops_GCT[[k]] = true.w
  print(paste("Done with cell = ", k, sep = ""))
}

names(cellpredictions_GCT) = paste("GCT_", cells, sep = "")
names(truecellprops_GCT) = paste("GCT_", cells, sep = "")

#save true cell proportionsm, predicted cell proportions, and LoD estimates for later
save(cellpredictions_GCT, truecellprops_GCT, LoDResults_127, file = "LoDResults_127.RData")


###############################################################################
# visualize/analyze the results
###############################################################################

par(mfrow = c(3,4))
for(k in 1:length(cells)) {
  LoB = LoBResults_127[1,cells[k]]
  LoD = LoDResults_127[1,cells[k]]
  tmp = hist(cellpredictions_GCT[[k]][, cells[k]]*100, xlab = "Cell Percentage (%)", 
             main = paste0("GCT_", cells[k]), xlim = c(0,3), cex.main = 2, cex.axis = 1.3, cex.lab = 1.5, freq = F, col = "grey90")
  abline(v = LoB, lty = "dotted", col = "red")
  abline(v = LoD, lty = "dotted", col = "blue")
  text(2.25 , .6, labels = paste("LoB = ", round(LoB,3), sep = ""), cex = 1, col = "red")
  text(2.25 , .6 +.25, labels = paste("LoD = ", round(LoD,3), sep = ""), cex = 1, col = "Blue") 
  
}

dev.off()







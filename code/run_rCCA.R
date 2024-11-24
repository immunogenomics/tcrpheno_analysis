
### This script is designed for mutliple parallel runs, each with a different
### parameter configuration

## That is, to launch 10 different parameter configurations:
#######for i in {1..10}; do
##########sbatch <SLURM flags> --wrap "Rscript run_rCCA.R $i"
#######done
  
############################ Parameter grid ###################################


# @l1s: lambda1 values to trial for ridge regularization of CCA
l1s = 10^(seq(-5,5))
# @l2s: lambda2 values to trial for ridge regularization of CCA
l2s = 10^(seq(-5,5))

# @perms: permute the rows of Y with respect to X? (used to assess statistical signifcance)
######## value of 1 = no permutation
######## value of 2+: permutation with a random seed
##e.g. perms = seq(1,11) will result in one rCCA without permutation, as well as 10 rCCAs following random permutation
perms = seq(1,11)


## all combinations of the above parameters are assembled into a grid,
## "Rscript run_rCCA.R $i" executes row $i of this grid
job_grid = expand.grid(l1s, l2s, perms)

args = commandArgs(trailingOnly = TRUE)
i = as.numeric(as.character(args[1]))
lambda1 = as.numeric(as.character(job_grid$Var1[i]))
lambda2 = as.numeric(as.character(job_grid$Var2[i]))
perm = as.numeric(as.character(job_grid$Var3[i]))


############################# other parameters ###########################

dir = ##<your output directory>

# @mode:
####### value "tuning": tuning lambda values though cross-validation within the training set
####### value "fitting": fitting the entire training dataset
mode = "tuning"

########################## helper functions ################################

trycatch_rcc <- function(X, Y, nCVs, lambda1, lambda2){
  fit = tryCatch(
    {
      mixOmics::rcc(X, Y, ncomp=nCVs, lambda1, lambda2, method="ridge")
    },
    error = function(cond) {
      message(cond)
    },
    finally={
      message("fit rcc")
    })
  return(fit)
}

get_permuted_order <- function(perm, barcodes, force_denovo=FALSE){
  samples = unique(train_sampleIDs)
  df = data.frame(cell = rownames(y), sample = train_sampleIDs)
  for (i in 1:length(samples)){
    cur = as.character(df$cell[df$sample==samples[i]])
    cur = sample(cur, length(cur), replace=FALSE)
    df$cell[df$sample==samples[i]] = cur
  }
  return(as.character(df$cell))
}

################################################################

library(Matrix)
library(dplyr)
library(mixOmics)
library(corpcor)
library(matrixcalc)

tag = paste0(c(mode, paste("lambda1.", lambda1, sep=""), paste("lambda2.", lambda2, sep="")), collapse="_")
if (perm>1){
  tag = paste(tag, paste("permutation", perm, sep=""), sep="_")
}
tag <- paste(tag, "rds", sep=".")

md = readRDS("data/merged_metadata_Datasets1and2.rds")

x = readRDS("data/CRtrtest_061324/CR_xtrain.rds")
y = readRDS("data/CRtrtest_061324/CR_ytrain.rds")

df = data.frame(cell = rownames(x))
df = left_join(df, md[,c("cell", "donor", "sample")])
train_sampleIDs <- as.character(df$donor)

mns_x = sapply(1:ncol(x), function(i) mean(x[,i], na.rm=TRUE))
names(mns_x) = colnames(x)
sds_x = sapply(1:ncol(x), function(i) sd(x[,i], na.rm=TRUE))
names(sds_x) = colnames(x)
mns_y = sapply(1:ncol(y), function(i) mean(y[,i], na.rm=TRUE))
names(mns_y) = colnames(y)
sds_y = sapply(1:ncol(y), function(i) sd(y[,i], na.rm=TRUE))
names(sds_y) = colnames(y)

for (j in 1:ncol(x)){
  x[,j] <- scale(x[,j])[,1]
}

for (j in 1:ncol(y)){
  y[,j] <- scale(y[,j])[,1]
}

ncomp = min(c(10, ncol(x), ncol(y)))
x[is.na(x)] <- 0

if (perm>1){
  rownames(y) = get_permuted_order(perm, rownames(x), force_denovo=TRUE)
  y = y[as.character(rownames(x)),]
}

setwd(dir)

if (mode=="tuning"){
  source("tune_rCCA_lambda.R")
  tuned = my.tune.rcc(as.matrix(x), as.matrix(y), grid1=c(lambda1), grid2=c(lambda2), folds=5, nCVs = 5, tag=tag, holdout="bypatient", metadata_file="data/merged_metadata_Datasets1and2.rds", patient_key="donor")
  saveRDS(tuned, file=paste("rCCA", tag, sep="_"))
} else if (mode=="fitting"){
  tag = paste(paste(lambda1, lambda2, sep="_"), tag, sep="_")
  res = NA
  res = trycatch_rcc(x, y, 10, lambda1, lambda2)
  while ((object.size(res)==object.size(NA)) & perm>1){
    rownames(y) = get_permuted_order(perm, rownames(x), force_denovo=TRUE)
    y = y[as.character(rownames(x)),]
    res = trycatch_rcc(x, y, 10, lambda1, lambda2)
  }
  res$X <- NULL
  res$Y <- NULL
  if (perm==1){
    save(res, file=paste("rCCA", tag, sep="_"))
  } else {
    saveRDS(res$cor, file=paste("rCCAcors", tag, sep="_"))
  }
} 
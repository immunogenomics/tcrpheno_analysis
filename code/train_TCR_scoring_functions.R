
scale_variables <- function(data, mns, sds){
  data = data[,colnames(data) %in% names(mns)]
  data = sweep(data, 2, mns)
  data = sweep(data, MARGIN=2, FUN="/", sds)
  data[is.na(data)] <- 0
  return(data)
}

library(doMC)
library(Matrix)
library(dplyr)
library(glmnet)

################################# Parameters #############################

# @cv: Which canonical variate to focus on?
######## value = 1: CV1, PLZF-high innate-like T cells
######## value = 2: CV2, CD8T vs CD4 T cells
######## value = 3: CV3, regulatory T cells
######## value = 4: CV4, memory T cells
cv = 1

# @chain: Which TCR chains to include in the scoring function?
######## value = "alphabeta": both TCRa and TCRb
######## value = "alpha": TCRa only
######## value = "beta": TCRb only
chain = "alphabeta"

tag = paste("target", cv, sep="")
tag = paste(tag, chain, sep="_")

M = readRDS("data/CRtrtest_061324/CR_xtrain.rds")

mns_x = sapply(1:ncol(M), function(x) mean(M[,x], na.rm=TRUE))
sds_x = sapply(1:ncol(M), function(x) sd(M[,x], na.rm=TRUE))
names(mns_x) = colnames(M)
names(sds_x) = colnames(M)
M = scale_variables(M, mns_x, sds_x)

if (chain=="beta"){
  M = M[,!(grepl("TRA", colnames(M)))]
  M = M[,!(grepl("TCRA", colnames(M)))]
}

if (chain=="alpha"){
  M = M[,!(grepl("TRB", colnames(M)))]
  M = M[,!(grepl("TCRB", colnames(M)))]
}

targets = read.csv("data/COMBATrenmerged_cellstatetargets_0606.csv")
targets = targets[targets$group=="train",]

if (cv==1){
  targets$target = targets$target1
}
if (cv==2){
  targets$target = targets$target2
}
if (cv==3){
  targets$target = targets$target3
}
if (cv==4){
  targets$target = targets$target4
}

M = M[rownames(M) %in% targets$cell[!(is.na(targets$target))],]
vars = sapply(1:ncol(M), function(x) var(M[,x]))
M = M[,vars!=0]

M = as.matrix(M)
rownames(targets) = as.character(targets$cell)
targets = targets[as.character(rownames(M)),]

pf = rep(1, ncol(M))

y = targets$target

registerDoMC(cores = 4)

fit = cv.glmnet(x=M, y=y, alpha=0, penalty.factor = pf, family="binomial", type.measure="mse", nfolds=5, parallel=TRUE)

saveRDS(fit, file=paste(tag, "LR_TCRscore.RData", sep="_"))



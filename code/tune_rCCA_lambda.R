
## modified source code from mixOmics::tune.rcc (https://rdrr.io/cran/mixOmics/man/tune.rcc.html)

####### modified for the hierarchical structure of single-cell data,
####### ensuring that cross-validation splits do not separate cells 
####### from the same individual

library(mixOmics)

trycatch_fit <- function(X, Y, omit, nCVs, lambda1, lambda2){
 fit = tryCatch(
                  {
                  rcc(X[-omit, , drop = FALSE], Y[-omit, , drop = FALSE], ncomp=nCVs, lambda1, lambda2, method="ridge")
                  },
                  error = function(cond) {
                      message(cond)
                     # return(NA)
                  },
                  finally={
                     message("fit fold")
                  })
  return(fit)
}

my.Mfold = function(X, Y, lambda1, lambda2, folds = 10, nCVs = 1, tag="defaulttag")
{
    posdef = TRUE
    print(paste("doing lambdas", paste(lambda1, lambda2)))   
    xscore = NULL
    yscore = NULL
    M = length(folds)
    for (m in 1:M)
    {
        print(paste("doing fold", m))
        omit = folds[[m]]
        X[omit, ][is.na(X[omit, ])] = 0
        Y[omit, ][is.na(Y[omit, ])] = 0
        fit = NA
        fit = trycatch_fit(X, Y, omit, nCVs, lambda1, lambda2) 
       if (object.size(fit)==object.size(NA)){ 
	  print("fit was NA")
          posdef = FALSE
          break
       }
       print("not NA fit") 
       nx = X[omit, , drop = FALSE] %*% fit$loadings$X[,1:nCVs]
       ny = Y[omit, , drop = FALSE] %*% fit$loadings$Y[,1:nCVs]
        xscore = rbind(xscore, nx)
        yscore = rbind(yscore, ny)
    }
    if (posdef){
      cv.score = sapply(1:ncol(xscore), function(x) cor(xscore[,x], yscore[,x]))
      print(cv.score)
      return(invisible(cv.score))
    } else {
      return(rep(NA, nCVs))
    }
}

my.tune.rcc <- function(X, Y, grid1, grid2, folds=10, nCVs = 1, tag="defaulttag", holdout="bycell", metadata_file="/raychaudhuri/data/ren/ren_metadata.rds", patient_key="PatientID", cell_key="cell"){

  validation = "Mfold"
  grid1 = unique(grid1)
  grid2 = unique(grid2)
  grid = expand.grid(grid1, grid2)
  nr = nrow(X)
  M = length(folds)
  if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > nr){
      stop("Invalid number of folds.")
  } else {
      M = round(folds)
      if (holdout=="bypatient") {
        df = data.frame(cell = as.character(rownames(Y)))
        md = readRDS(metadata_file)
        md$cell = as.character(md[,which(colnames(md)==cell_key)])
        md$PatientID = as.character(md[,which(colnames(md)==patient_key)])
        df = left_join(df, md[,c("cell", "PatientID")])
        mtx_pt = as.character(df$PatientID)
        pts = unique(mtx_pt)
        sfolds = split(sample(pts), rep(1:M, length=length(pts)))
        folds = list()
        for (m in 1:M){
           folds[[m]] = which(mtx_pt %in% sfolds[[m]])
           names(folds)[m] = m
         }
      } else {
        folds = split(sample(1:nr), rep(1:M, length = nr))
      }
  }
  cv.score = matrix(nrow=nrow(grid), ncol=nCVs)
  for (i in 1:nrow(grid)){
     cv.score = rbind(cv.score, my.Mfold(X=X, Y=Y, lambda1 = as.numeric(as.character(grid$Var1[i])), lambda2 = as.numeric(as.character(grid$Var2[i])), folds=folds, nCVs = nCVs, tag=tag)) 
  }
   cv.score.grid = cbind(grid, cv.score)
   cv.score.grid = cv.score.grid[complete.cases(cv.score.grid),]
   cv.score = sapply(1:nrow(cv.score.grid), function(x) sum(cv.score.grid[x,3:ncol(cv.score.grid)]))
   opt = cv.score.grid[cv.score == max(cv.score, na.rm=TRUE),1:2]
   if (nrow(opt)==0){
      print(cv.score.grid)
      print("no matrices invertible, all NAs")
      return(NA)
   }
   out = list(opt.lambda1 = opt[[1]], opt.lambda2 = opt[[2]], 
                 opt.score = max(cv.score.grid), grid1 = grid1, grid2 = grid2, score.grid = cv.score.grid, scores = cv.score)  
   return(out)
}


## Function Description:
##   Calculate Cluster Robust Standard Errors from Fitted Model Objects
##
## TO-DO ITEMS: 1. generalize to other regression objects: e.g., "lm", "zeroinfl"
##              2. calculate a sandwich estimate of theta (overdispersion term)

clrobustit <- function(fit.model, dv, clusterid, small.sample=FALSE) {
  require(lmtest, quietly=TRUE)
  require(sandwich, quietly=TRUE)
  
  ## fitted model class
  model.class <- class(fit.model)
  
  sndwch <- function(model,family,rank,y,id,vcv,estfn,small.n=small.sample) {
    ## degree of freedom adjustment robust variance estimator
    dfc.adj <- function(M,N,K,dist=family, smalln=small.n) {
      dfc0 <- (M/(M-1)) * ((N-1)/(N-K))  # Gaussian (Stata)
      dfc1 <- M/(M-1)                    # Non-Gaussian (Stata)
      dfc2 <- M/(M-K)                    # Mancl & DeRouen
      
      if (smalln) {
        return(dfc2)
      } else if (dist=="gaussian") {
        return(dfc0)
      } else {
        return(dfc1)
      }
    }
    
    N.obs <- length(y)
    dfc <- dfc.adj(M=length(unique(id)), N=N.obs, K=length(rank))
    
    if (missing(vcv)) vcv <- vcov(model)
    if (missing(estfn)) estfn <- estfun(model)
    
    uj <- apply(estfn, 2, function(x) tapply(x, id, sum))
    
    bread <- N.obs * vcv
    meat  <- crossprod(uj)/N.obs
    
    vcovCL <- dfc*(1/N.obs * (bread %*% meat %*% bread))
    
    return(vcovCL)
  }
  
  if (identical(model.class,"lm") | identical(model.class,c("glm","lm"))) {
    ## ascertain the distribution type
    dist <- ifelse(identical(model.class,"lm"), "gaussian", fit.model$family$family)
    
    vcovCL <- sndwch(model=fit.model, family=dist, rank=fit.model$rank, y=dv,
                     id=clusterid, small.n=small.sample)
    
    ## perform robust Wald tests in each submodel
    res <- coeftest(fit.model, vcov=vcovCL)

  } else if (model.class=="hurdle" | model.class=="zeroinfl") {
    require(pscl, quietly=TRUE)
    
    ## ascertain the family of each submodel
    family.bin <- fit.model$dist$zero
    family.cnt <- fit.model$dist$count

    ## create separate binary and truncated count variables
    dv.rnd <- round(dv)
    dv.bin <- as.numeric(ifelse(dv.rnd > 0, 1, 0))
    dv.cnt <- as.numeric(ifelse(dv.rnd > 0, dv.rnd, NA))
    
    ## create temporary data frames to determine cluster sizes for hurdle submodels
    dv.df  <- data.frame(cbind(clusterid, dv.cnt, dv.bin))
    dv.bin.df <- dv.df[!is.na(dv.bin), ]
    dv.cnt.df <- dv.df[!is.na(dv.cnt), ]
    
    ## total parameters (K) in each submodel
    K.bin <- length(fit.model$optim$zero$par)     # binary model
    K.cnt <- length(fit.model$optim$count$par)    # truncated count model
        
    ## Assemble bread and meat for sandwich estimator
    vcov.all <- fit.model$vcov
    vcov.bin <- vcov.all[grep("^zero_",  rownames(vcov.all), value = FALSE), grep("^zero_",  colnames(vcov.all), value = FALSE)]
    vcov.cnt <- vcov.all[grep("^count_", rownames(vcov.all), value = FALSE), grep("^count_", colnames(vcov.all), value = FALSE)]
    
    estfun.all <- estfun(fit.model)
    estfun.bin <- estfun.all[!is.na(dv.df$dv.bin), grep("^zero_",  colnames(estfun.all), value = FALSE)]
    estfun.cnt <- estfun.all[!is.na(dv.df$dv.cnt), grep("^count_", colnames(estfun.all), value = FALSE)]
    
    vcovCL.bin <- sndwch(model=fit.model, family=family.bin, rank=K.bin, y=dv.bin.df$dv.bin,
                         id=dv.bin.df$clusterid, vcv=vcov.bin, estfn=estfun.bin)
    vcovCL.cnt <- sndwch(model=fit.model, family=family.bin, rank=K.cnt, y=dv.cnt.df$dv.cnt,
                         id=dv.cnt.df$clusterid, vcv=vcov.cnt, estfn=estfun.cnt)

    # Perform robust Wald tests in each submodel
    res <- list(bin = coeftest(fit.model, vcov=vcovCL.bin),
                cnt = coeftest(fit.model, vcov=vcovCL.cnt)
                )
  } else {
    stop("The provided model object is not currently supported.")  
  }

  # Output robust result(s)
  return(res)
}

### Description:
##   Calculates cluster robust standard error estimates from fixed-effect GLMs.
##   This function supports R's built-in glm() and lm() functions as well as
##   hurdle and zero-inflated count models from the pscl package.
##
###  Author: David Huh
##
###  Dependencies: lmtest, sandwich
##
###  Arguments:  fit.model  = a fitted model object of class glm, lm, hurdle, or zeroinfl
##               clusterid  = a vector defining which cluster each observation belongs to.
##                            Per Miglioretti & Heagerty (2007), in fully-nested designs,
##                            clustering should be defined at the top level.
##             small.sample = apply the small-sample degrees-of-freedom correction suggested
##                            by Mancl & DeRouen (2001)
##
###  Values:       sims = a matrix of simulations (rows: simulates, columns: predictor combo)
##              summary = a regression summary based on cluster robust standard errors
##                 coef = a vector of the mean coeficient estimates
##                 vcov = the cluster robust variance-covariance matrix
##
### References:
##
##    1. Mancl, L. A., & DeRouen, T. A. (2001). A Covariance estimator for GEE with improved
##         small-sample properties. Biometrics, 57, 126–134.
##
##    2. Miglioretti, D. L., & Heagerty, P. J. (2007). Marginal modeling of nonnested multilevel
##         data using standard software. American Journal of Epidemiology, 165, 453–463.
##         doi: 10.1093/aje/kwk020
##

clrobustse <- function(fit.model, clusterid, small.sample=FALSE) {
  require(lmtest, quietly=TRUE)
  require(sandwich, quietly=TRUE)
   
  model.class <- class(fit.model)  # fitted model class
  
  ## generic function for calculating clustered sandwich estimates
  sndwch <- function(model,family,rank,id,vcv,estfn,small.n) {
    
    ## degree of freedom adjustment robust variance estimator
    dfc.adj <- function(M,N,K,dist,smalln) {
      dfc0 <- (M/(M-1)) * ((N-1)/(N-K))  # Gaussian (Stata)
      dfc1 <- M/(M-1)                    # Non-Gaussian (Stata)
      dfc2 <- M/(M-K)                    # Mancl & DeRouen
      
      if (smalln) return(dfc2) else
        if (dist=="gaussian") return(dfc0) else 
          return(dfc1)
      }
    }
    
    N.obs <- length(id)            # number of observations
    N.clust <- length(unique(id))  # number of clusters
    
    dfc <- dfc.adj(M=N.clust, N=N.obs, K=rank, dist=family, smalln=small.n)
    
    if (missing(vcv)) vcv <- vcov(model)
    if (missing(estfn)) estfn <- estfun(model)
    
    uj <- apply(estfn, 2, function(x) tapply(x, id, sum))
    
    bread <- N.obs * vcv
    meat  <- crossprod(uj)/N.obs
    
    vcovCL <- dfc*(1/N.obs * (bread %*% meat %*% bread))
    
    return(vcovCL)
  }
  
  if (identical(model.class, "lm")) {
    dist <- "gaussian"
    
    vcovCL <- sndwch(model=fit.model, family=dist, rank=fit.model$rank,
                     id=clusterid, small.n=small.sample)
    
    ## bundle result for 'lm' model
    res <- list(summary=coeftest(fit.model, vcov=vcovCL),
                coef=coef(fit.model),
                vcov=vcovCL)
    
  } else if (identical(model.class, c("glm","lm"))) {
    ## ascertain the distribution type
    dist <- fit.model$family$family
    
    vcovCL <- sndwch(model=fit.model, family=dist, rank=fit.model$rank,
                     id=clusterid, small.n=small.sample)

    ## bundle result for 'glm' model
    res <- list(summary=coeftest(fit.model, vcov=vcovCL),
                coef=coef(fit.model),
                vcov=vcovCL)
    
  } else if (model.class=="hurdle" | model.class=="zeroinfl") {
    require(pscl, quietly=TRUE)
    
    ## ascertain the family of each submodel
    family.bin <- fit.model$dist$zero
    family.cnt <- fit.model$dist$count

    ## create separate binary and truncated count variables
    dv.rnd <- round(fit.model$y)
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
    
    vcovCL.bin <- sndwch(model=fit.model, family=family.bin, rank=K.bin,
                         id=dv.bin.df$clusterid, vcv=vcov.bin, estfn=estfun.bin, small.n=small.sample)
    vcovCL.cnt <- sndwch(model=fit.model, family=family.bin, rank=K.cnt,
                         id=dv.cnt.df$clusterid, vcv=vcov.cnt, estfn=estfun.cnt, small.n=small.sample)

    ## Extract estimated regression coefficients
    coef.all <- coef(fit.model)
    coef.bin <- coef.all[grep("^zero_",  names(coef.all), value = FALSE)]
    coef.cnt <- coef.all[grep("^count_", names(coef.all), value = FALSE)]
    
    ## Bundle submodel output
    res.bin <- list(summary=coeftest(fit.model, vcov=vcovCL.bin),
                    coef=coef.bin,
                    vcov=vcovCL.bin)
    
    res.cnt <- list(summary=coeftest(fit.model, vcov=vcovCL.cnt),
                    coef=coef.cnt,
                    vcov=vcovCL.cnt)
    
    res <- list(bin.model = res.bin,
                cnt.model = res.cnt)
    
  } else {
    stop("The model object provided is not currently supported.")  
  }

  # Output result(s)
  return(res)
}

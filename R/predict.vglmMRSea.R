predict.vglmMRSea <- function(object, newdata=NULL, newdists=NULL, type="response", coeff=NULL, includeB0=TRUE, conf_int=FALSE) {

  library(mvtnorm)
  # If newdata and newdists are included predict from them
  # Otherwise predict from object
  if (is.null(newdata)) {
    data <- object@data
  } else {
    data <- newdata
  }
  
  splineParams <- object@splineParams
  
  # Create bases for newdata - 2d smooth
  if (is.null(newdists)){
    datDist <- splineParams[[1]]$dist
  } else {
    datDist <- newdists
  }
  
  aR <- splineParams[[1]]$knotPos
  radii <- splineParams[[1]]$radii
  radiusIndices <- splineParams[[1]]$radiusIndices
  #varz <- object@varshortnames
  varz <- names(object@misc$orig.assign)
  varz_short <- object@varshortnames
  n_sp <- length(splineParams)
  
  if (varz[1]=="(Intercept)"){
    intercept_present <- TRUE
    if (length(varz) > 1){
      varz <- varz[2:length(varz)]
    } else {
      varz <- c()
    }
  } else {
    intercept_present <- FALSE
  }
  
  if (length(varz) > 0){
    if (substr(varz[length(varz)],1,5)=="LRF.g"){
      Twod_present <- TRUE
      if (length(varz) > 1){
        varz <- varz[1:(length(varz)-1)]
      } else {
        varz <- c()
      }
    } else {
      Twod_present <- FALSE
    }
  } else {
    Twod_present <- FALSE
  }
  
  if (length(varz) > 0) {
    # Create bases for newdata - 1d variables first
    for (nv in 1:length(varz)) {
      var_in_txt <- paste0(varz_short[nv], "<-data$", varz_short[nv])
      var_in <- eval(parse(text=var_in_txt))
      #var_bs_txt <- paste0("bs(", varz_short[nv], ", knots=splineParams[[", nv + 1,
      #                     "]]$knots, degree=splineParams[[", nv + 1, "]]$degree, ",
      #                     "Boundary.knots=splineParams[[", nv+ 1, "]]$bd)")
      var_bs_txt <- varz[nv]
      var_bs <- eval(parse(text=var_bs_txt))
      if (nv==1){
        new_bs_all <- var_bs
      } else {
        new_bs_all <- cbind(new_bs_all, var_bs)
      }
    }
  }

  bs_lrf <- LRF.g(radiusIndices, datDist, radii, aR)
  
  if (Twod_present){
    if (length(varz) > 0) {
      new_bs_all <- cbind(new_bs_all, bs_lrf)
    } else {
      new_bs_all <- bs_lrf
    }
  }
  
  if(is.null(coeff)) {
    coefs <- object@coefficients
  } else {
    coefs <- coeff
  }
  
  if (intercept_present) {
    if (exists("new_bs_all")){
      new_bs_all <- cbind(rep(1, nrow(new_bs_all)), new_bs_all)
    } else {
      new_bs_all <- rep(1, nrow(datDist))
    }
  }
  
  #nyy <- dim(object@y)[2] 
  #neta <- nyy - 1
  neta <- dim(object@predictors)[2]
  
  coefs_mat <- matrix(coefs, ncol=neta, byrow=TRUE)
  
  preds_out = new_bs_all %*% coefs_mat
  
  if (type == "response") {
    preds_out = object@family@linkinv(preds_out, extra=object@extra)
  }
  
  if (conf_int == TRUE) {
    sum_out <- summaryvglm(object)
    covs <- sum_out@cov.unscaled
    rcoefs <- rmvnorm(1000, coefs, covs)
    #quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}
    #cis <- apply(rcoefs, 2, quant.func)
    #lw_ci <- matrix(cis[1,], ncol=neta, byrow=T)
    #hi_ci <- matrix(cis[2,], ncol=neta, byrow=T)
    #lw_lim <- new_bs_all %*% lw_ci
    #hi_lim <- new_bs_all %*% hi_ci
    #if (type == "response") {
    #  lw_lim = object@family@linkinv(lw_lim, extra=object@extra)
    #  hi_lim = object@family@linkinv(hi_lim, extra=object@extra)
    #}
    get_bootstrap <- function(rcoef){
      bsin <- matrix(rcoef, ncol=neta, byrow=T)
      linkvals <- new_bs_all %*% bsin
      respvals <- object@family@linkinv(linkvals, extra=object@extra)
      return(respvals)
    }
    resvalz <- get_bootstrap(rcoefs[1,])
    # try storing as separate arrays for each outcome more efficient to add to smaller arrays
    yvalz <- rep(list(matrix(NA, nrow=nrow(resvalz), ncol=1000)),ncol(resvalz))
    for (ii in 1:1000){
      resval <- get_bootstrap(rcoefs[ii,])
      for (yy in 1:ncol(resvalz)){
        yvalz[[yy]][,ii] <- resval[,yy]
      }
    }
    lowerci <- apply(yvalz[[1]], 1, quantile, probs=c(0.025))
    upperci <- apply(yvalz[[1]], 1, quantile, probs=c(0.975))
    
    for (aa in 2:6) {
      yval <- yvalz[[aa]]
      anicilw <- apply(yval, 1, quantile, probs=c(0.025))
      anicihi <- apply(yval, 1, quantile, probs=c(0.975))
      lowerci <- cbind(lowerci, anicilw)
      upperci <- cbind(upperci, anicihi)
    }
    namez <- colnames(object@y)
    colnames(lowerci) <- namez
    colnames(upperci) <- namez
    return(list("predictions"=preds_out, "lower_limit"=lowerci, "higher_limit"=upperci))
    
  } else {
    return(preds_out)
  }

}



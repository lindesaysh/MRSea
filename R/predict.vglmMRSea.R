predict.vglmMRSea <- function(object, newdata=NULL, type="response", coeff=NULL, newdists=NULL, includeB0=TRUE) {
  
  if (is.null(newdata)) {
    data <- object@data
  } else {
    data <- newdata
  }
  
  splineParams <- object@splineParams
  aR <- splineParams[[1]]$knotPos
  radii <- splineParams[[1]]$radii
  radiusIndices <- splineParams[[1]]$radiusIndices
  varz <- object@varshortnames
  
  if (length(varz) > 0) {
    # Create bases for newdata - variables first
    for (nv in 1:length(varz)) {
      var_in_txt <- paste0(varz[nv], "<-data$", varz[nv])
      var_in <- eval(parse(text=var_in_txt))
      var_bs_txt <- paste0("bs(", varz[nv], ", knots=splineParams[[", nv+1,
                           "]]$knots, degree=splineParams[[", nv+1, "]]$degree, ",
                           "Boundary.knots=splineParams[[", nv+1, "]]$bd)")
      var_bs <- eval(parse(text=var_bs_txt))
      if (nv==1){
        new_bs_all <- var_bs
      } else {
        new_bs_all <- cbind(new_bs_all, var_bs)
      }
    }
  }
  
  # Create bases for newdata - 2d smooth
  if (is.null(newdists)) {
    distMats <- makeDists(
      cbind(data$x.pos, data$y.pos),
      na.omit(splineParams[[1]]$knotgrid)
    )
    datDist <- distMats$dataDist
  } else {
    datDist <- newdists
  }
  bs_lrf <- LRF.g(radiusIndices, datDist, radii, aR)
  
  if (length(varz) > 0) {
    new_bs_all <- cbind(new_bs_all, bs_lrf)
  } else {
    new_bs_all <- bs_lrf
  }
  
  if(is.null(coeff)) {
    coefs <- object@coefficients
  } else {
    coefs <- coeff
  }
  
  if (includeB0) {
    new_bs_all <- cbind(rep(1, nrow(new_bs_all)), new_bs_all)
  }
  
  nyy <- dim(object@y)[2] 
  neta <- nyy - 1
  
  coefs <- matrix(coefs, ncol=neta, byrow=TRUE)
  
  preds_out = new_bs_all %*% coefs
  
  if (type == "response") {
    preds_out = object@family@linkinv(preds_out, extra=object@extra)
  }
  
  return(preds_out)
}



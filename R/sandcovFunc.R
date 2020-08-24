sandcov<-function (model, id=NULL){
  if(is.null(id)){
    id<-1:length(model$y)
  }
  model$weights <- pmax(1e-08, model$weights)
  l <- length(coef(model))
  vv <- summary.glm(model)$cov.unscaled
  mm <- model.matrix(model)
  ii <- match(colnames(vv), colnames(mm))
  infl <- (residuals(model, "response") * model$prior.weights * 
             mm[, ii]) %*% vv
  jj <- match(colnames(mm), colnames(vv))
  infl <- infl[, jj]
  colnames(infl) <- colnames(mm)
  infl <- rowsum(infl, id)
  t(infl) %*% infl
}

clsandcov   <- function(dat,fm, cluster){
  #]attach(dat, warn.conflicts = F)
  library(sandwich)
  
  if(min(table(cluster))==0){
    cluster<-droplevels(cluster)
    warning("Length zero found for 1 or more clusters, droplevels() has been applied")
  }
  
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  return(vcovCL)}

sandcovS4<-function (model, id=NULL){
  if(is.null(id)){
    id<-1:length(model@y)
  }
  model@weights <- pmax(1e-08, model@weights)
  l <- length(coef(model))
  vv <- summary.glm(model)@cov.unscaled
  mm <- model.matrix(model)
  ii <- match(colnames(vv), colnames(mm))
  infl <- (residuals(model, "response") * model@prior.weights * 
             mm[, ii]) %*% vv
  jj <- match(colnames(mm), colnames(vv))
  infl <- infl[, jj]
  colnames(infl) <- colnames(mm)
  infl <- rowsum(infl, id)
  t(infl) %*% infl
}


clsandcov_mn   <- function(dat,fm, cluster){
 # attach(dat, warn.conflicts = F)
  library(sandwich)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm@rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  return(vcovCL)
}


bread.gamMRSea<-function (x, ...) 
{
  if (!is.null(x$na.action)) 
    class(x$na.action) <- "omit"
  sx <- summary.glm(x)
  wres <- as.vector(residuals(x, "working")) * weights(x, 
                                                       "working")
  dispersion <- if (substr(x$family$family, 1, 17) %in% c("poisson", 
                                                          "binomial", "Negative Binomial")) 
    1
  else sum(wres^2)/sum(weights(x, "working"))
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])) * dispersion)
}
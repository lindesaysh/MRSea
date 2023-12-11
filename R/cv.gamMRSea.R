#' Cross-validation for gamMRSea Models
#' 
#' This function calculates the estimated K-fold cross-validation prediction error for gamMRSea models. It has been adapted from the \link[boot]{cv.glm} function from the \code{boot} library
#' 
#' @param data A matrix or data frame containing the data. The rows should be cases and the columns correspond to variables, one of which is the response.
#' @param modelobject An object of class "gamMRSea" containing the results of a gamMRSea model fitted to data.
#' @param cost A function of two vector arguments specifying the cost function for the cross-validation. The first argument to cost should correspond to the observed responses and the second argument should correspond to the predicted or fitted responses from the generalized linear model. cost must return a non-negative scalar value. The default is the average squared error function.
#' @param K The number of groups into which the data should be split to estimate the cross-validation prediction error. The value of K must be such that all groups are of approximately equal size. If the supplied value of K does not satisfy this criterion then it will be set to the closest integer which does and a warning is generated specifying the value of K used. The default is to set K equal to the number of observations in data which gives the usual leave-one-out cross-validation.
#' @param replicate (\code{default=FALSE}).  When using the replicate function, and panels are specified in the model object, \code{replicate=TRUE} will change the seed and select new panel based folds for each iteration.
#' @param s.eed (\code{default = NULL}). If \code{NULL} then a seed is randomly generated and stored in the attributes of the output.  If specified, the value for \code{s.eed} is used and also stored in the attributes of the output. 
#' 
#' @details For more information please see the \link[boot]{cv.glm} function in the boot library
#' 
#' @examples
#' 
#' # load data
#' data(ns.data.re)
#' ns.data.re$foldid<-getCVids(ns.data.re, folds=5)
#'  
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#' cv.gamMRSea(data=ns.data.re, modelobject = model, K=5)$delta
#' 
#' @export
#' 
#'

cv.gamMRSea<-function (data, modelobject, cost = function(y, yhat) mean((y - yhat)^2), 
          K = n, replicate=FALSE, s.eed=NULL) 
{
  
  if(!is.null(modelobject$panels)){
    if(length(unique(modelobject$panels)) < K){
      stop("Not enough unique panels to make K folds. Please reduce K and try again")
    }
  }
    
  
  call <- match.call()
  # if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
  #   s.eed<-sample(1:100000, size = 1)
  #   set.seed(s.eed)
  # seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  # 
  if(is.null(s.eed)){
    s.eed<-sample(1:100000, size = 1)
  }
  set.seed(s.eed)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  n <- nrow(data)
  if ((K > n) || (K <= 1)) 
    stop("'K' outside allowable range")
  K.o <- K
  K <- round(K)
  kvals <- unique(round(n/(1L:floor(n/2))))
  temp <- abs(kvals - K)
  if (!any(temp == 0)) 
    K <- kvals[temp == min(temp)][1L]
  if (K != K.o) 
    warning(gettextf("'K' has been set to %f", K), domain = NA)
  f <- ceiling(n/K)
  
  if(replicate==TRUE){
      if(length(unique(modelobject$panels))==nrow(data) | is.null(modelobject$panels)){
        # have indep
        s <- boot:::sample0(rep(1L:K, f), n)  
      }else{
        # not indep
        #s.eed<-sample(1:100000, size = 1)
        s<-getCVids(data, K, block=modelobject$panels, seed = s.eed)  
      }
  }else{
    if(is.null(modelobject$cvfolds)){
      s <- boot:::sample0(rep(1L:K, f), n)  
    }else{
      if(length(table(modelobject$cvfolds)) == K){
        s<-modelobject$cvfolds
      }else{
        #s.eed<-sample(1:100000, size = 1)
        s<-getCVids(data, K, block=modelobject$panels, seed = s.eed)
        warning("K in CV call is not the same as K used for cvfolds stored in model object. \n New folds have been created using s.eed parameter (see attributes of output) \n Note: the model object folds have not been changed")
      }
    }  
  }
  
  
  n.s <- table(s)
  glm.y <- modelobject$y
  cost.0 <- cost(glm.y, fitted(modelobject))
  ms <- max(s)
  CV <- 0
  Call <- modelobject$call
  for (i in seq_len(ms)) {
    
    j.out <- seq_len(n)[(s == i)]
    j.in <- seq_len(n)[(s != i)]
    
    
    Call$data <- data[j.in, , drop = FALSE]
    
    if(!is.null(modelobject$splineParams)){
      splineParams<-modelobject$splineParams
      if(!is.null(splineParams[[1]]$dist)){
        splineParams[[1]]$dist<-modelobject$splineParams[[1]]$dist[j.in,]
        g2k<-modelobject$splineParams[[1]]$dist[j.out,]
        Call$splineParams<-splineParams
      }else{
        g2k<-NULL
        Call$splineParams<-splineParams
      }
    }else{
      splineParams<-NULL
      g2k<-NULL
    }
    
    
    d.glm <- eval.parent(Call)
    p.alpha <- n.s[i]/n
    cost.i <- cost(glm.y[j.out], predict(object=d.glm, newdata=data[j.out, ,drop = FALSE], g2k=g2k, type = "response"))
    CV <- CV + p.alpha * cost.i
    
    if(!is.null(modelobject$splineParams[[1]]$dist)){
      g2k<-modelobject$splineParams[[1]]$dist
    }else{
      g2k<-NULL
    }
    
    cost.0 <- cost.0 - p.alpha * cost(glm.y, predict(object=d.glm, 
                                                     newdata=data, g2k=g2k, type = "response"))
  }
  
  CV0 <- CV + cost.0
  
  if(exists("CV0")){
    if(CV0 < 0){
      CV0 <- Inf}
  }
  
  delta <- as.numeric(c(CV, CV0))
  
  output <- list(call = call, 
                 K = K,
                 delta = delta)
  
  attr(output, "s.eed") <- s.eed
  attr(output, "seed") <- seed
  
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cv.gamMRSea.mn<-function (data, modelobject, cost = function(y, yhat) mean((y - yhat)^2), K = n, replicate=FALSE) {
  call <- match.call()
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  n <- nrow(data)
  if ((K > n) || (K <= 1)) {
    stop("'K' outside allowable range")
  }
  K.o <- K
  K <- round(K)
  kvals <- unique(round(n/(1L:floor(n/2))))
  temp <- abs(kvals - K)
  if (!any(temp == 0)) {
    K <- kvals[temp == min(temp)][1L]
  }
  if (K != K.o) {
    warning(gettextf("'K' has been set to %f", K), domain = NA)
  }
  f <- ceiling(n/K)
  
  ### function only used for multinomials so model object must be S4 check if is vglm or vglmMRSea
  ### if not vglmMRSea can't contain panels etc already so must be indep
  if(replicate==TRUE){
    if (class(modelobject)[1]=="vglm") {
      s <- boot:::sample0(rep(1L:K, f), n)  
    } else {
      if(length(unique(modelobject@panels))==nrow(data) | is.null(modelobject@panels)){
        # have indep
        s <- boot:::sample0(rep(1L:K, f), n)  
      }else{
        # not indep
        s.eed<-sample(1:100000, size = 1)
        ### I don't think there is anything in getcvids that means it won't work with S4
        s<-getCVids(data, K, block=modelobject@panels, seed = s.eed) 
      }
    }
  }else{
    if (class(modelobject)[1]=="vglm") {
      s <- boot:::sample0(rep(1L:K, f), n)  
    } else {
      if(length(modelobject@cvfolds)<1){
        s <- boot:::sample0(rep(1L:K, f), n)  
      }else{
        s<-modelobject@cvfolds
      } 
    }
  }
  ### being lazy for now just check if cvfolds exists or not if it is still base vglm class it can't exist
  if (class(modelobject)[1]=="vglm"){
    splineParams <- NULL
  } else {
    splineParams<-modelobject@splineParams
  }
  
  n.s <- table(s)
  glm.y <- modelobject@y
  cost.0 <- cost(glm.y, fitted(modelobject))
  ms <- max(s)
  CV <- 0
  Call <- modelobject@call
  for (i in seq_len(ms)) {
    
    j.out <- seq_len(n)[(s == i)]
    j.in <- seq_len(n)[(s != i)]
    
    
    Call$data <- data[j.in, , drop = FALSE]
    
    if(class(modelobject)[1]=="vglmMRSea"){
      if (length(modelobject@splineParams)>0) {
        splineParams <- modelobject@splineParams
        if(!is.null(splineParams[[1]]$dist)){
          splineParams[[1]]$dist<-splineParams[[1]]$dist[j.in,]
          g2k<-splineParams[[1]]$dist[j.out,]
          Call$splineParams<-splineParams
        }else{
          g2k<-NULL
          Call$splineParams<-splineParams
        }
      } else {
        splineParams<-NULL
        g2k<-NULL
      }
    }else{
      splineParams<-NULL
      g2k<-NULL
    }
    
    
    d.glm <- eval.parent(Call)
    p.alpha <- n.s[i]/n
    cost.i <- cost(glm.y[j.out,], predict(object=d.glm, newdata=data[j.out, ,drop = FALSE], g2k=g2k, type = "response"))
    CV <- CV + p.alpha * cost.i
    
    if(!is.null(splineParams[[1]]$dist)){
      g2k<-splineParams[[1]]$dist
    }else{
      g2k<-NULL
    }
    
    cost.0 <- cost.0 - p.alpha * cost(glm.y, predict(object=d.glm, 
                                                     newdata=data, g2k=g2k, type = "response"))
  }
  list(call = call, K = K, delta = as.numeric(c(CV, CV + cost.0)), 
       seed = seed)
}
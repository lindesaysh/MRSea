#' Cross-validation for gamMRSea Models
#' 
#' This function calculates the estimated K-fold cross-validation prediction error for gamMRSea models. It has been adapted from the \link[boot]{cv.glm} function from the \code{boot} library
#' 
#' @param data A matrix or data frame containing the data. The rows should be cases and the columns correspond to variables, one of which is the response.
#' @param modelobject An object of class "gamMRSea" containing the results of a gamMRSea model fitted to data.
#' @param cost A function of two vector arguments specifying the cost function for the cross-validation. The first argument to cost should correspond to the observed responses and the second argument should correspond to the predicted or fitted responses from the generalized linear model. cost must return a non-negative scalar value. The default is the average squared error function.
#' @param K The number of groups into which the data should be split to estimate the cross-validation prediction error. The value of K must be such that all groups are of approximately equal size. If the supplied value of K does not satisfy this criterion then it will be set to the closest integer which does and a warning is generated specifying the value of K used. The default is to set K equal to the number of observations in data which gives the usual leave-one-out cross-validation.
#' @param replicate (\code{default=FALSE}).  When using the replicate function, and panels are specified in the model object, \code{replicate=TRUE} will change the seed and select new panel based folds for each iteration.
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
          K = n, replicate=FALSE) 
{
  call <- match.call()
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
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
        s.eed<-sample(1:100000, size = 1)
        s<-getCVids(data, K, block=modelobject$panels, seed = s.eed)  
      }
  }else{
    if(is.null(modelobject$cvfolds)){
      s <- boot:::sample0(rep(1L:K, f), n)  
    }else{
      s<-modelobject$cvfolds
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
  list(call = call, K = K, delta = as.numeric(c(CV, CV + cost.0)), 
       seed = seed)
}
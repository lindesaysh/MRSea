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

cv.vglmMRSea<-function (modelobject, datain=NULL, distsin=NULL, rawdistsin=NULL, K = n, cost=cost_mn_accuracy, replicate=FALSE) 
{
  
  splineParams <- modelobject@splineParams
  if (is.null(distsin)){
    datDist <- splineParams[[1]]$dist
  } else {
    datDist <- distsin
  }
  if (is.null(rawdistsin)){
    raw_dists <- modelobject@splineParams[[1]]$raw_dists
  } else {
    raw_dists <- rawdistsin
  }
  if (is.null(datain)) {
    data <- modelobject@data
    glm.y <- modelobject@y
    fitted.y <- fitted(modelobject)
  } else {
    data <- datain
    y_as_factor <- datain$response
    levs <- levels(y_as_factor)
    for (lv in 1:length(levs)) {
      level_mask <- y_as_factor == levs[lv]
      if (lv == 1) {
        glm.y = level_mask
      } else {
        glm.y = cbind(glm.y, level_mask)
      }
    }
    colnames(glm.y) <- levs
    fitted.y <- predict.vglmMRSea(object=d.glm, newdata=data, newdists=datDist, type = "response")
  }
  
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
    if(length(unique(modelobject@panels))==nrow(data) | is.null(modelobject@panels)){
      # have indep
      s <- boot:::sample0(rep(1L:K, f), n)  
    }else{
      # not indep
      s.eed<-sample(1:100000, size = 1)
      s<-getCVids(data, K, block=modelobject@panels, seed = s.eed)  
    }
  }else{
    if(length(modelobject@cvfolds)==0){
      s <- boot:::sample0(rep(1L:K, f), n)  
    }else{
      s<-modelobject@cvfolds
    }
  }
  n.s <- table(s)
  ms <- max(s)
  Call <- modelobject@call
  
  CV <- 0
  cost.0 <- cost_mn_accuracy(glm.y, fitted.y)
  
  for (i in seq_len(ms)) {
    
    j.out <- seq_len(n)[(s == i)]
    j.in <- seq_len(n)[(s != i)]
    
    Call$data <- data[j.in, , drop = FALSE]
    
    splineParams<-modelobject@splineParams
    splineParams[[1]]$dist<-modelobject@splineParams[[1]]$dist[j.in,]
    Call$splineParams<-splineParams

    d.glm <- eval(Call)
    
    p.alpha <- n.s[i]/n
    predz <- predict.vglmMRSea(object=d.glm, newdata=data[j.out,], newdists=datDist[j.out,], type = "response")
    resps <- glm.y[j.out,]
    cost.i <- cost_mn_accuracy(resps, predz)
    CV <- CV + p.alpha * cost.i
    
    predz2 <- predict.vglmMRSea(object=d.glm, newdata=data[j.in,], newdists=datDist[j.in,], type = "response")
    resps2 <- glm.y[j.in,]
    cost.0 <- cost.0 - p.alpha * cost_mn_accuracy(resps2, predz2)
  }

  list(call = call, K = K, delta = as.numeric(c(CV, CV + cost.0)), 
       seed = seed)
}


cost_mn_accuracy <- function(response, predictions){
  resp_class <- apply(response,1,which.max)
  pred_class <- apply(predictions,1,which.max)
  correct_preds <- pred_class == resp_class
  accuracy <- sum(correct_preds) / length(correct_preds)
  cost <- 1- accuracy
}
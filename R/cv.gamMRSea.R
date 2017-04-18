


cv.gamMRSea<-function (data, glmfit, cost = function(y, yhat) mean((y - yhat)^2), 
          K = n) 
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
  s <- boot:::sample0(rep(1L:K, f), n)
  n.s <- table(s)
  glm.y <- glmfit$y
  cost.0 <- cost(glm.y, fitted(glmfit))
  ms <- max(s)
  CV <- 0
  Call <- glmfit$call
  for (i in seq_len(ms)) {
    
    j.out <- seq_len(n)[(s == i)]
    j.in <- seq_len(n)[(s != i)]
    
    
    Call$data <- data[j.in, , drop = FALSE]
    
    if(!is.null(glmfit$splineParams)){
      splineParams<-glmfit$splineParams
      if(!is.null(splineParams[[1]]$dist)){
        splineParams[[1]]$dist<-glmfit$splineParams[[1]]$dist[j.in,]
        g2k<-glmfit$splineParams[[1]]$dist[j.out,]
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
    
    if(!is.null(glmfit$splineParams[[1]]$dist)){
      g2k<-glmfit$splineParams[[1]]$dist
    }else{
      g2k<-NULL
    }
    
    cost.0 <- cost.0 - p.alpha * cost(glm.y, predict(object=d.glm, 
                                                     newdata=data, g2k=g2k, type = "response"))
  }
  list(call = call, K = K, delta = as.numeric(c(CV, CV + cost.0)), 
       seed = seed)
}
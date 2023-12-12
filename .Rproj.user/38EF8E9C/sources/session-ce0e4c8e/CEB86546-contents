

cumres.gamMRSea<-function (object, variable = c("predicted", colnames(
  model.matrix(object))), 
          data = data.frame(model.matrix(object)), R = 1000, b = 0, 
          plots = min(R, 50), breakties = 1e-12, seed = round(runif(1, 
                                                                    1, 1e+09)), ...) 
{
  dginv <- function(z) 1
  a.phi <- 1
  switch(class(object)[1], lm = {
    a.phi <- summary(object)$sigma^2
    h <- function(z) 1
  }, glm = {
    f <- family(object)
    if (f$family == "gaussian") {
      a.phi <- summary(object)$dispersion
    }
    g <- f$linkfun
    ginv <- f$linkinv
    dginv <- f$mu.eta
    canonf <- do.call(f$family, list())
    caninvlink <- canonf$linkinv
    canlink <- canonf$linkfun
    Dcaninvlink <- canonf$mu.eta
    Dcanlink <- function(x) 1/Dcaninvlink(canlink(x))
    h <- function(z) Dcanlink(ginv(z)) * dginv(z)
    covmat<-summary(object)$cov.scaled
  }, gamMRSea = {
    f <- family(object)
    if (f$family == "gaussian") {
      a.phi <- summary(object)$dispersion
    }
    g <- f$linkfun
    ginv <- f$linkinv
    dginv <- f$mu.eta
    canonf <- do.call(f$family, list())
    caninvlink <- canonf$linkinv
    canlink <- canonf$linkfun
    Dcaninvlink <- canonf$mu.eta
    Dcanlink <- function(x) 1/Dcaninvlink(canlink(x))
    h <- function(z) Dcanlink(ginv(z)) * dginv(z)
    covmat<-summary(object)$cov.robust
  },stop("Unsupported model!"))
  response <- all.vars(formula(object))[1]
  X <- model.matrix(object)
  n <- nrow(X)
  r <- residuals(object, type = "response")
  yhat <- predict(object=object, type = "response")
  beta <- coef(object)
  if (any(is.na(beta))) 
    stop("Over-parametrized model")
  Xbeta <- X %*% beta
  etaraw <- (as.numeric(dginv(Xbeta)) * X)
  hatW.MC <- function(x) {
    myorder <- order(x)
    x <- x[myorder]
    Ii <- covmat
    A <- as.vector(h(Xbeta) * r)/a.phi
    S <- apply(X, 2, function(x) x * A)
    beta.iid <- Ii %*% t(S[myorder, , drop = FALSE])
    r0 <- r[myorder]
    D0 <- etaraw[myorder, , drop = FALSE]
    Wfun <- "W2"
    if (b != 0) {
      Wfun <- "W"
    }
    output <- .C(Wfun, R = as.integer(R), b = as.double(b), 
                 n = as.integer(n), npar = as.integer(nrow(Ii)), 
                 xdata = as.double(x), rdata = as.double(r0), betaiiddata = as.double(beta.iid), 
                 etarawdata = as.double(D0), plotnum = as.integer(plots), 
                 seed = as.integer(seed), KS = as.double(0), CvM = as.double(0), 
                 Wsd = as.double(numeric(n)), cvalues = as.double(numeric(R)), 
                 Ws = as.double(numeric(plots * n)), Wobs = as.double(numeric(n)), 
                 pkg = "gof")
    return(list(output = output, x = x))
  }
  if (!is.na(match(response, variable))) 
    variable[match(response, variable)] <- "predicted"
  variable <- unique(variable)
  UsedData <- X[, na.omit(match(variable, colnames(X))), drop = FALSE]
  myvars <- colnames(UsedData)[apply(UsedData, 2, function(x) length(unique(x)) > 
                                       2)]
  if ("predicted" %in% variable) 
    myvars <- c("predicted", myvars)
  untie <- runif(n, 0, breakties)
  W <- c()
  What <- c()
  Wsd <- c()
  KS <- c()
  CvM <- c()
  mytype <- c()
  UsedVars <- c()
  UsedData <- c()
  allcvalues <- c()
  for (v in myvars) {
    x <- NULL
    if (v == "predicted") {
      x <- yhat
    }
    else if (v %in% colnames(X)) {
      x <- X[, v]
    }
    if (!is.null(x)) {
      UsedVars <- c(UsedVars, v)
      onesim <- hatW.MC(x + untie)
      UsedData <- cbind(UsedData, onesim$x)
      W <- cbind(W, onesim$output$Wobs)
      Wsd <- cbind(Wsd, onesim$output$Wsd)
      What <- c(What, list(matrix(onesim$output$Ws, nrow = n)))
      KS <- c(KS, onesim$output$KS)
      CvM <- c(CvM, onesim$output$CvM)
      allcvalues <- cbind(allcvalues, onesim$output$cvalues)
      mytype <- c(mytype, "residual")
    }
    else cat("Variable '", v, "' not found.\n", sep = "")
  }
  if (length(UsedVars) < 1) 
    return(NULL)
  res <- list(W = W, What = What, x = UsedData, KS = KS, CvM = CvM, 
              R = R, n = nrow(UsedData), sd = Wsd, cvalues = allcvalues, 
              variable = UsedVars, type = mytype, untie = untie, model = class(object)[1])
  class(res) <- "cumres"
  res
}
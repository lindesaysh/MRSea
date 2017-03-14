#' gamMRSea model function
#' 
#' See \code{\link{glm}} for details.  A \code{splineParams} object may be specified as part of the model object.
#' 
#' @author LAS Scott-Hayward, University of St Andrews
#' 
#' @export
#' 
gamMRSea<-function (formula, family = gaussian, data, weights, subset,
                    na.action, start = NULL, etastart, mustart, offset, control = list(...),
                    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, splineParams=NULL, ...)
{
  
  if(!is.null(splineParams)){
    dists<-splineParams[[1]]$dist
    aR<-splineParams[[1]]$invInd[splineParams[[1]]$knotPos]
    radii<-splineParams[[1]]$radii
    radiusIndices<-splineParams[[1]]$radiusIndices
  
    # create new dataframe:
    
    }

  call <- match.call()
  
    if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset", "splineParams"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(MRSea::model.frame.gamMRSea)
  mf <- eval(mf, envir=parent.frame())
  if (identical(method, "model.frame.gamMRSea"))
    return(mf)
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start, etastart = etastart,
                   mustart = mustart, offset = offset, family = family,
                   control = control, intercept = attr(mt, "intercept") >
                     0L))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
                      x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights,
                      offset = offset, family = family, control = control,
                      intercept = TRUE))
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  
  varshortnames<-names(data)
  
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
                     data = data, offset = offset, control = control, method = method, contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf), splineParams=splineParams, varshortnames=varshortnames))
  class(fit) <- c(fit$class, c("gamMRSea", "glm", "lm"))
  fit
}

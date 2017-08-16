#' Summarising model fits from models fitted using the MRSea package.
#'
#' @param object	an object of class "gamMRSea", usually, a result of a call from the MRSea package.
#' @param x	an object of class "summary.gamMRSea", usually, a result of a call to summary.gamMRSea.
#' @param dispersion	the dispersion parameter for the family used. Either a single numerical value or NULL (the default), when it is inferred from object (see 'Details').
#' @param correlation	logical; if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param digits	the number of significant digits to use when printing.
#' @param symbolic.cor	logical. If TRUE, print the correlations in a symbolic form (see symnum) rather than as numbers.
#' @param signif.stars	logical. If TRUE, 'significance stars' are printed for each coefficient.
#' @param varshortnames vector stating the short versions of the covariate names if required.
#' @param ...	further arguments passed to or from other methods.
#'
#' @details
#' \code{print.summary.gamMRSea} tries to be smart about formatting the coefficients, standard errors, etc. and additionally gives 'significance stars' if signif.stars is TRUE. The coefficients component of the result gives the estimated coefficients and their estimated standard errors (raw and robust), together with their ratio (from robust s.e.). The third column gives the robust standard errors calculated using the sandwich estimator.  If no correlation is present, the second and third columns are the same as the sandwich estimator is not used when data points are independent. The fourth column is labelled Wald and gives the Wald test statistic, based on the robust standard errors. The fifth column gives the two-tailed p-value corresponding to the Wald test ().
#'
#' Aliased coefficients are omitted in the returned object but restored by the print method.
#'
#' Correlations are printed to two decimal places (or symbolically): to see the actual correlations print summary(object)$correlation directly.
#'
#' \code{summary.gamMRSea} returns an object of class "summary.gamMRSea", a list with components
#' call	the component from object.
#' family	the component from object.
#' deviance	the component from object.
#' contrasts the component from object.
#' df.residual	the component from object.
#' null.deviance	the component from object.
#' df.null	the component from object.
#' deviance.resid	the deviance residuals: see residuals.glm.
#' coefficients	the matrix of coefficients, standard errors, z-values and p-values. Aliased coefficients are omitted.
#' aliased	named logical vector showing if the original coefficients are aliased.
#' dispersion	either the supplied argument or the inferred/estimated dispersion if the latter is NULL.
#' df	a 3-vector of the rank of the model and the number of residual degrees of freedom, plus number of coefficients (including aliased ones).
#' cov.unscaled	the unscaled (dispersion = 1) estimated covariance matrix of the estimated coefficients.
#' cov.scaled	ditto, scaled by dispersion.
#' correlation	(only if correlation is true.) The estimated correlations of the estimated coefficients.
#' symbolic.cor	(only if correlation is true.) The value of the argument symbolic.cor.
#'
#' @author Lindesay Scott-Hayward, Univeristy of St Andrews.
#' @note Code adapted from \code{summary.glm}
#'
#' @examples 
#' 
#' # load data
#' data(ns.data.re)
#' ns.data.re$foldid<-getCVids(ns.data.re, folds=5)
#'  
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#' summary(model)
#'
#' @export

summary.gamMRSea<-function (object, dispersion = NULL, varshortnames=NULL, ...)
{
  if(is.null(object$panels)){panelid<-1:nrow(object$data)
  }else{
    panelid<-object$panels
  }

  if(!is.null(object$varshortnames)){
    varshortnames=object$varshortnames
    object<-summaryshortnames(object, varshortnames)
  }

  vbeta<-sandcov(object, panelid)

  

  est.disp <- FALSE
  df.r <- object$df.residual
  if (is.null(dispersion))
    dispersion <- if (object$family$family %in% c("poisson",
                                                  "binomial"))
      1
  else if (df.r > 0) {
    est.disp <- TRUE
    if (any(object$weights == 0))
      warning("observations with zero weight not used for calculating dispersion")
    sum((object$weights * object$residuals^2)[object$weights >
                                                0])/df.r
  }
  else {
    est.disp <- TRUE
    NaN
  }
  aliased <- is.na(coef(object))
  p <- object$rank
  if (p > 0) {
    p1 <- 1L:p
    Qr <- qr.lm(object)
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
    covmat <- dispersion * covmat.unscaled
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    
    if(max(table(panelid))==1){
      s.err.sand <-sqrt(var.cf)
      vbeta<-covmat
    }else{
      s.err.sand <- sqrt(diag(vbeta))  
    }
    
    tvalue <- coef.p/s.err.sand
    traw <- coef.p/s.err
    dn <- c("Estimate", "Std. Error", "Robust S.E.")
    if (!est.disp) {
      pvalue <- 2 * pnorm(-abs(tvalue))
      rawp <- 2 * pnorm(-abs(traw))
      coef.table <- cbind(coef.p, s.err, s.err.sand,tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                    "z value", "Pr(>|z|)"))
    }
    else if (df.r > 0) {
      pvalue <- 2 * pt(-abs(tvalue), df.r)
      rawp <- 2 * pt(-abs(traw), df.r)
      coef.table <- cbind(coef.p, s.err,  s.err.sand, tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                    "t value", "Pr(>|t|)"))
    }
    else {
      coef.table <- cbind(coef.p, NaN, NaN, NaN, NaN)
      dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                    "t value", "Pr(>|t|)"))
    }
    df.f <- NCOL(Qr$qr)
  }
  else {
    coef.table <- matrix(, 0L, 4L)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                         "Robust Std. Err.",
                                         "t value", "Pr(>|t|)"))
    covmat.unscaled <- covmat <- matrix(, 0L, 0L)
    df.f <- length(aliased)
  }
  keep <- match(c("call", "terms", "family", "deviance", "aic",
                  "contrasts", "df.residual", "null.deviance", "df.null",
                  "iter", "na.action"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object,
                                                         type = "deviance"), coefficients = coef.table, aliased = aliased,
                              dispersion = dispersion, df = c(object$rank, df.r, df.f),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat, cov.robust=vbeta, panelid=panelid, rawpvals=rawp))
  class(ans) <- "summary.gamMRSea"
  return(ans)
}


print.gamMRSea<-function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts))
      cat("  [contrasts: ", apply(cbind(names(co), co),
                                  1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
      x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Null Deviance:\t   ", format(signif(x$null.deviance,
                                           digits)), "\nResidual Deviance:", format(signif(x$deviance,
                                                                                           digits)), "\tAIC:", format(signif(x$aic, digits)))
  cat("\n")
  invisible(x)
}



print.summary.gamMRSea<-function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
          signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Deviance Residuals: \n")
  if (x$df.residual > 5) {
    x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                          na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
  }

  xx <- zapsmall(x$deviance.resid, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df <- if ("df" %in% names(x))
      x[["df"]]
    else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
          sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                               colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
      format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
                                                                "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
                                                                                                                                 "deviance")]), digits = max(5L, digits + 1L)), " on",
                                                 format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                                           1L, paste, collapse = " "), sep = "")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)))

  if(max(table(x$panelid))==1){
    maxpanels <- paste(max(table(x$panelid)), ' (independence assumed)', sep='')
  }else{maxpanels<-max(table(x$panelid))}
  
  cat("\n\nMax Panel Size = ", maxpanels,"; Number of panels = ", length(unique(x$panelid)),"\nNumber of Fisher Scoring iterations: ", x$iter, "\n", sep = "")


  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2L), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}


qr.lm<-function (x, ...)
{
  if (is.null(r <- x$qr))
    stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  r
}


#-----------------------------------------------------------------------------
#' Calculate percentile confidence intervals from a matrix of bootstrapped predictions
#' 
#' @param preds matrix of bootstrap predictions where each column is a bootstrap realisation
#' @param quants (\code{default = c(0.025, 0.975)}. Vector of length two of quantiles. 
#' 
#' @examples 
#' \dontrun{
#' makeBootCIs(bootPreds)
#' }
#' @export
#' 

makeBootCIs<-function(preds, quants=c(0.025,0.975)){
  cis<-t(apply(preds, 1,  quantile, probs = quants, na.rm=T ))
}
  
  